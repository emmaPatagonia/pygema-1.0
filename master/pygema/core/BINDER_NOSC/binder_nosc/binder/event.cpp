/*
 * This file is part of binder.
 *
 * Copyright (C) 2009 Andy Heath, Stuart Nippress & Andreas Rietbrock,
 *                    University of Liverpool
 *
 * This work was funded as part of the NERIES (JRA5) project.
 * Additional funding for Nippress from NERC research grant NE/C000315/1
 *
 * binder is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * binder is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with binder.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <iomanip>

#include "utl_time.h"
#include "event.h"
#include "debug.h"

#include "livloc.h"
#include <hypo71.h>
#include <aehnll.h>

using namespace Andy;

//debug
static void printPicks(std::list<MyPick> picks)
{
  for (std::list<MyPick>::iterator iter = picks.begin(); iter != picks.end(); iter++) {
    std::cerr << *iter << " [" << iter->relocatedWeight << "]" << std::endl;
  }
	std::cerr << "----------" << std::endl;
}

static bool checkEventWindow(std::list<MyPick>& pickList, double pickTime, float eventWindow)
{
  if (!pickList.empty()) {
    MyPick earliestPick = pickList.front();
    if (fabs(earliestPick.pickTime - pickTime) > eventWindow) {
      return false;
    }
  }

  return true;  
}

//--------------------------------------------------

static bool checkForRelatedP(const MyPick& pick, std::list<MyPick>& pickList, MyPick& pPick)
{
  std::list<MyPick>::iterator miter;
  std::string stationCode = pick.stationCode;
  
  for (miter = pickList.begin(); miter != pickList.end(); miter++) {
    if (miter->stationCode == stationCode) {
      if (miter->phaseHint == std::string("P")) {
        pPick = *miter;
        return true;
      }
    }
  }
  
  return false;
}

//--------------------------------------------------

Event::Event(double _originTime, const Position& _location, std::list<MyPick> _pickList)
: originTime(_originTime), location(_location), rms(-9999.0f), gap(-9999.0f), numberP(0), numberS(0), pickList(_pickList)
{
  pickList.sort();
}

//--------------------------------------------------

Event::Event(double _originTime, float _rms, const Position& _location, std::list<MyPick> _pickList)
: originTime(_originTime), location(_location), rms(_rms), gap(-9999.0f), numberP(0), numberS(0), pickList(_pickList)
{
  pickList.sort();
}

//--------------------------------------------------

// order by time
bool Event::operator <(const Event& event)
{
  return originTime < event.originTime;
}

//--------------------------------------------------

bool Event::operator ==(const Event& event)
{
  return originTime == event.originTime && location == event.location;
}

//--------------------------------------------------

// add a single pick
void Event::addPick(const MyPick& pick)
{
  pickList.push_back(pick);
  pickList.sort();
}

//--------------------------------------------------

// if the pick was generated by this event, store the pick
bool Event::associatePick(const Position& stationPosition,
                          float (*ttFunc)(const Position&, const Position&, float&),
                          float pCutoff, float sCutoff,
                          float pEventWindow, float sEventWindow,
                          float pToSVelocityRatio,
                          MyPick pick)
{
  double pt, et, tt, pRes, sRes;
  bool possibleP = false;
  bool possibleS = false;
  MyPick dummy;
	float ainDummy;

  pt = pick.pickTime;

  // if we have any picks yet, check time between earliest pick and new one is within sEventWindow
  if (!checkEventWindow(pickList, pt, sEventWindow)) return false;

  // travel time of P wave from this event to station
  tt = ttFunc(location, stationPosition, ainDummy);
  
  // check for P
  et = originTime + tt;
  if ((pRes = fabs(pt - et)) <= pCutoff) {
    possibleP = true;
  }
  
  // check for S
  if (checkForRelatedP(pick, pickList, dummy)) { // this here because hypo71 needs one P for each S
    et = originTime + (pToSVelocityRatio * tt);
    if ((sRes = fabs(pt - et)) <= sCutoff) {
      possibleS = true;
    }
  }
  
  if (possibleS) {
    pick.phaseHint = std::string("S");
  }
  else if (possibleP) {
    pick.phaseHint = std::string("P");
  }
  else {
    return false;
  }
  
  //pick.historyList.push_back(1); // history
  
  pick.relocatedWeight = 1111.0f; //debug

  addPick(pick);

  return true;
}

//--------------------------------------------------

// remove the picks in this event that are also in 'primary'... if there are enough picks
// left over (based on eventThreshold) then this is still an event. If there are less, this
// event will be allowed to die and the left over picks are now unassociated
bool Event::filterPicks(unsigned int eventThreshold, Event& primary, std::list<MyPick>& unassociatedPickList)
{
  std::list<MyPick> remainingPickList;
  std::list<MyPick>::iterator thisIter, primaryIter;
  bool foundFlag;

  for (thisIter = pickList.begin(); thisIter != pickList.end(); thisIter++) {
    foundFlag = false;
    for (primaryIter = primary.pickList.begin(); primaryIter != primary.pickList.end(); primaryIter++) {
      if (*thisIter == *primaryIter) {
        foundFlag = true;
        break;
      }
    }
    if (!foundFlag) {
      remainingPickList.push_back(*thisIter);
    }
  }

  if (remainingPickList.size() >= eventThreshold) {
    pickList = remainingPickList;
    return true;
  }
  else {
    // temp
    /*
    for (thisIter = remainingPickList.begin(); thisIter != remainingPickList.end(); thisIter++) {
      thisIter->historyList.push_back(3); // history
    }
    */
    // end temp
    unassociatedPickList.merge(remainingPickList);
    return false;
  }
}

/*
 * DEBUG
 *
static void debugHypoPickList(std::list<MyPick>& pickList, std::list<Hypo71::HypoPick>& hypoPickList)
{
  std::list<Hypo71::HypoPick>::iterator hiter;
  std::list<MyPick>::iterator miter;
  std::string stationCode1, phase1, stationCode2, phase2;
  
  if (hypoPickList.empty()) {
    std::cerr << "EMPTY!" << std::endl;
  }
  else {
    for (hiter = hypoPickList.begin(); hiter != hypoPickList.end(); hiter++) {
      stationCode2.clear();
      phase2.clear();
      stationCode1 = hiter->pickPtr->waveformID().stationCode();
      phase1 = hiter->pickPtr->phaseHint();
			int refCnt = hiter->pickPtr->referenceCount();
      if (hiter->pPickPtr) {
        stationCode2 = hiter->pPickPtr->waveformID().stationCode();
        phase2 = hiter->pPickPtr->phaseHint();
      }
      std::cerr << stationCode1 << " " << phase1 << " " << refCnt << " " << stationCode2 << " " << phase2 << std::endl;
    }
  }
  
  std::cerr << std::endl << "Picklist size = " << pickList.size() << std::endl;
  for (miter = pickList.begin(); miter != pickList.end(); miter++) {
    phase1 = miter->pickPtr->phaseHint();
    stationCode1 = miter->pickPtr->waveformID().stationCode();
    std::cerr << stationCode1 << " " << phase1 << std::endl;
  }
  std::cerr << "------------------------------" << std::endl;
}

static bool debugRepeatPicks(std::list<MyPick> pickList)
{
  std::string testCode, testPhase, tryCode, tryPhase;
  std::list<MyPick>::iterator miter;
  
  for (miter = pickList.begin(); miter != pickList.end(); miter++) {
    tryCode = miter->pickPtr->waveformID().stationCode();
    tryPhase = miter->pickPtr->phaseHint();
    if (testCode == tryCode && testPhase == tryPhase) {
      return true;
    }
    else {
      testCode = tryCode;
      testPhase = tryPhase;
    }
  }
  
  return false;
}
*
* END DEBUG
*/

//--------------------------------------------------

bool Event::relocateHypo71(const Cartesian& cart)
{
  Hypo71 hypo71;
  Hypo71::Sol sol;
  std::list<Hypo71::HypoPick> hypoPickList;
  std::list<Hypo71::HypoPick>::iterator hiter;
  std::list<MyPick>::iterator miter;
  MyPick pPick;

  // create list of pick objects used in relocation
  for (miter = pickList.begin(); miter != pickList.end(); miter++) {
    Hypo71::HypoPick hypoPick;
    hypoPick.stationCodeStr = miter->stationCode;
    hypoPick.pickTime = miter->pickTime;
    hypoPick.phaseHintStr = miter->phaseHint;
    hypoPick.importedWeight = miter->weight;
    if (miter->polarity == POLARITY_POSITIVE) {
      hypoPick.polarityStr = std::string("U");
    }
    else if (miter->polarity == POLARITY_NEGATIVE) {
      hypoPick.polarityStr = std::string("D");
    }
    else {
      hypoPick.polarityStr = std::string("_");
    }

    if (miter->phaseHint == std::string("P")) {
      hypoPickList.push_back(hypoPick);
    }
    else {
      if (checkForRelatedP(*miter, pickList, pPick)) {
        hypoPick.hasRelatedPPick = true;
        hypoPickList.push_back(hypoPick);
      }
    }
  }

  if (hypo71.locate(hypoPickList, sol)) {
    originTime = sol.orig_time;
    cart.toXYZ(sol.event_lon, sol.event_lat, &(location.x), &(location.y));
    location.z = sol.event_depth;
    error = Position(sol.lat_error, sol.lon_error, sol.depth_error);
    rms = sol.rms_residual;
    gap = sol.gap;
    numberP = sol.num_of_p;
    numberS = sol.num_of_s;
    // copy the residuals from relocation objects into pickList
    for (hiter = hypoPickList.begin(), miter = pickList.begin(); hiter != hypoPickList.end(); hiter++, miter++) {
      miter->timeResidual = hiter->residual;
      miter->relocatedWeight = hiter->weight;
    }
  }
  else {
    std::cerr << "Event::relocate: call to hypo71 failed" << std::endl;
    return false;
  }
  
  return true;
}

//--------------------------------------------------
/*
bool Event::relocateLivLoc( float largeResidual, const Cartesian& cart )
{
	std::list<MyPick>::iterator miter;
	std::list<LivLoc::LLPick> livLocPickList;
	std::list<LivLoc::LLPick>::iterator llIter;
	LivLoc::LLOrigin origin;

  // create list of pick objects used in relocation
  for (miter = pickList.begin(); miter != pickList.end(); miter++) {
		LivLoc::LLPick::Polarity polarity;
		if (miter->polarity == POLARITY_POSITIVE) {
      polarity = LivLoc::LLPick::POSITIVE;
    }
    else if (miter->polarity == POLARITY_NEGATIVE) {
      polarity = LivLoc::LLPick::NEGATIVE;
    }
    else {
      polarity = LivLoc::LLPick::UNDECIDABLE;
    }
		LivLoc::LLPick llPick = LivLoc::LLPick(miter->stationCode, miter->pickTime, miter->phaseHint, polarity, miter->weight);
		llPick.index = miter->index; // hack
    livLocPickList.push_back(llPick);
		// set the residual for the pick larger that acceptable in Event::relocate()
		miter->timeResidual = largeResidual;
	}

  livLocObj->relocate( livLocPickList, origin, false );

	originTime = origin.originTime;
	cart.toXYZ(origin.location.y, origin.location.x, &(location.x), &(location.y));
	location.z = origin.location.z;
	error = origin.error;
	rms = origin.rms;
	gap = origin.gap;
	numberP = numberS = 0;

	// we are only interested in keeping the 'arrivals' in this event
	for (llIter = origin.pickList.begin(); llIter != origin.pickList.end(); llIter++) {
		if ( llIter->isArrival ) {
			for (miter = pickList.begin(); miter != pickList.end(); miter++) {
				if ( miter->index == llIter->index ) { // relates to same pick in event pickList
					miter->timeResidual = llIter->timeResidual;
					miter->relocatedWeight = llIter->relocatedWeight;
					miter->distance = llIter->distance;
					miter->azimuth = llIter->azimuth;
          if ( miter->isSPhase() ) {
						numberS++;
					}
					else {
						numberP++;
					}
					break;
				}
			}
		}
	}
	
	return true;
}
*/
//--------------------------------------------------

bool Event::relocateLivLoc( const Cartesian& cart, std::list<MyPick>& unassociatedPickList )
{
	std::list<MyPick>::iterator miter;
	std::list<LivLoc::LLPick> livLocPickList;

  // create list of pick objects used in relocation
  for (miter = pickList.begin(); miter != pickList.end(); miter++) {
		LivLoc::LLPick::Polarity polarity;
		if (miter->polarity == POLARITY_POSITIVE) {
      polarity = LivLoc::LLPick::POSITIVE;
    }
    else if (miter->polarity == POLARITY_NEGATIVE) {
      polarity = LivLoc::LLPick::NEGATIVE;
    }
    else {
      polarity = LivLoc::LLPick::UNDECIDABLE;
    }
		LivLoc::LLPick llPick = LivLoc::LLPick(miter->stationCode, miter->pickTime, miter->phaseHint, polarity, miter->weight);
		llPick.index = miter->index; // hack
    livLocPickList.push_back(llPick);
	}

	std::list<LivLoc::LLPick>::iterator llIter;
	LivLoc::LLOrigin origin;

	// suggest an initial origin (need to concert back to lat/lon first)
	float lat, lon;
	cart.toLonLat(location.x, location.y, &lon, &lat);
  //origin.location = Position(lat, lon, -location.z);
  origin.location = Position(lat, lon, location.z);
  origin.originTime = originTime;

  if ( livLocObj->relocate( livLocPickList, origin, true ) ) { // the 'true' tells livLoc to start at 'origin.location'
    originTime = origin.originTime;
    cart.toXYZ(origin.location.y, origin.location.x, &(location.x), &(location.y));
		// location.z = -origin.location.z;
		location.z = origin.location.z;
		error = origin.error;
    rms = origin.rms;
    gap = origin.gap;
    numberP = numberS = 0;
		pickList.clear();
		// we are only interested in keeping the 'arrivals' in this event
		for (llIter = origin.pickList.begin(); llIter != origin.pickList.end(); llIter++) {
			Polarity polarity = POLARITY_UNDECIDABLE;
			if (llIter->polarity == LivLoc::LLPick::POSITIVE) {
				polarity = POLARITY_POSITIVE;
			}
			else if (llIter->polarity == LivLoc::LLPick::NEGATIVE) {
				polarity = POLARITY_NEGATIVE;
			}
			else {
				polarity = POLARITY_UNDECIDABLE;
			}
			MyPick pick(llIter->stationCode, llIter->pickTime, llIter->phaseHint, "Z", polarity, llIter->weight);
			pick.phaseHint = llIter->phaseHint; // added 24/6/10
			pick.timeResidual = llIter->timeResidual;
			pick.relocatedWeight = llIter->relocatedWeight;
			pick.distance = llIter->distance;
			pick.azimuth = llIter->azimuth;
			pick.index = llIter->index; // hack
			pick.originTime = originTime;
			if ( llIter->isArrival ) {
				pickList.push_back(pick);
				if ( pick.isSPhase() ) {
					numberS++;
				}
				else {
					numberP++;
				}
			}
			else {
				unassociatedPickList.push_back(pick);
			}
		}
  }
  else {
    //std::cerr << "Event::relocate: call to LivLoc::relocate failed" << std::endl;
    return false;
  }

	return true;
}

//--------------------------------------------------

static bool lessThanStationNames(const MyPick& mp1, const MyPick& mp2)
{
  return mp1.stationCode < mp2.stationCode;
}

//--------------------------------------------------

static bool compareStationNames(const MyPick& mp1, const MyPick& mp2)
{
  return mp1.stationCode == mp2.stationCode;
}

//--------------------------------------------------

// check and remove similar phases from the same station - we assume all phases are P here
bool Event::removeSimilarPicks()
{
  unsigned int oldSize = pickList.size();

  pickList.sort(lessThanStationNames);
  pickList.unique(compareStationNames);

  return oldSize != pickList.size();
}

//--------------------------------------------------

// relocate this event using a location algorithm (if this returns false it means the event no longer contains enough picks)
bool Event::relocate(const Cartesian& cart, float pResidualCutoff, float sResidualCutoff, unsigned int eventThreshold, std::list<MyPick>& unassociatedPickList)
{
  // new: to ensure the P's come before their respective S's
  pickList.sort();

	/*
	float largeResidual = pResidualCutoff > sResidualCutoff ?  pResidualCutoff : sResidualCutoff;
  if ( !relocateLivLoc( 2.0f * largeResidual, cart ) ) {
		*/
	if ( !relocateLivLoc( cart, unassociatedPickList ) ) {
  /*
    for (miter = pickList.begin(); miter != pickList.end(); miter++) {
      unassociatedPickList.push_back(*miter);
    }
   */
    //return false;
    //std::cerr << "relocateLivLoc returned false" << std::endl; // debug
  }

	/*
  if (!relocateHypo71(cart)) {
    return false;
  }
	*/
  /*
  if (!relocateNLL(cart, "nllpickfile.txt")) {
    return false;
  }
  */
  
  //printPicks(pickList); // debug

  std::list<MyPick> tempPickList;
  std::list<MyPick>::iterator miter;
  float residualCutoff;
  int tmp_no_picks;

	// remove any picks with poor residuals
  for (miter = pickList.begin(); miter != pickList.end(); miter++) {
    if (miter->phaseHint == std::string("P")) {
      residualCutoff = pResidualCutoff;
    }
    else {
      residualCutoff = sResidualCutoff;
    }
    if (fabsf(miter->timeResidual) < residualCutoff) {
      tempPickList.push_back(*miter);
    }
    else {
      //miter->historyList.push_back(-1); // history
      unassociatedPickList.push_back(*miter);
    }
  }

  // check if this event is still sustainable
  pickList = tempPickList;
  // now check if there are still enough observations close by e.g. less than 20s
  // now check if there are still enough observations close by e.g. less than 60s
  tmp_no_picks=0;
  for (miter = pickList.begin(); miter != pickList.end(); miter++) {
	if((miter->pickTime-miter->originTime) < 60.0) {
	    tmp_no_picks++;
	}
  }
  //if (pickList.size() < eventThreshold) {
  if (tmp_no_picks < eventThreshold) {
    // too few picks - add the remaining list contents to the unassociated ones
    for (miter = pickList.begin(); miter != pickList.end(); miter++) {
      //miter->historyList.push_back(-2); // history
      unassociatedPickList.push_back(*miter);
    }
    return false;
  }
  
  return true;;
}

//--------------------------------------------------

// non-class

static void updatePicksWithNLL(std::list<AehNLL::NLLPick>& nllPickList, Event *eptr)
{
  std::list<AehNLL::NLLPick>::iterator niter;
  std::list<MyPick>::iterator miter;

  eptr->numberP = 0;
  eptr->numberS = 0;

  for (niter = nllPickList.begin(), miter = eptr->pickList.begin(); niter != nllPickList.end(); niter++, miter++) {
    if (niter->phaseHintStr == std::string("P")) {
      eptr->numberP++;
    }
    else if (niter->phaseHintStr == std::string("S")) {
      eptr->numberS++;
    }
    miter->timeResidual = niter->residual;
    miter->relocatedWeight = niter->weight;
  }
}

//--------------------------------------------------

static void updateLocationWithNLL(const Cartesian& cart, const AehNLL::Sol& sol, Event *eptr)
{
  TIME utl;
  
  utl.yr = sol.year;
  utl.mo = sol.month;
  utl.day = sol.day;
  utl.hr = sol.hour;
  utl.mn = sol.min;
  utl.sec = float(sol.sec);

  eptr->originTime = base_diff(utl);
 
  cart.toXYZ(sol.dlong, sol.dlat, &(eptr->location.x), &(eptr->location.y));
  eptr->location.z = sol.depth;
  // incorrect mapping if ellipsod is oblique to axes
  eptr->error = Position(sol.ellipsoid[AehNLL::Sol::LEN1], sol.ellipsoid[AehNLL::Sol::LEN2], sol.ellipsoid[AehNLL::Sol::LEN3]);
  eptr->rms = sol.rms;
  eptr->gap = sol.gap;
  // number of P & S set in 'updatePicksWithNLL'
}

//--------------------------------------------------

bool Event::relocateNLL(const Cartesian& cart, const char *pickFilename)
{
  AehNLL aehNLL;
  std::list<AehNLL::NLLPick> nllPickList;
  AehNLL::Sol sol;
  std::list<MyPick>::iterator miter;
  
  for (miter = pickList.begin(); miter != pickList.end(); miter++) {
    AehNLL::NLLPick nllPick;
    nllPick.stationCodeStr = miter->stationCode;
    nllPick.pickTime = miter->pickTime;
    nllPick.phaseHintStr = miter->phaseHint;
    nllPick.methodIDStr = miter->methodID;
    if (miter->polarity == POLARITY_POSITIVE) {
      nllPick.polarityStr = std::string("U");
    }
    else if (miter->polarity == POLARITY_NEGATIVE) {
      nllPick.polarityStr = std::string("D");
    }
    else {
      nllPick.polarityStr = std::string("?");
    }
    nllPickList.push_back(nllPick);
  }
  
  if (!aehNLL.locate(pickFilename, nllPickList, sol)) {
    std::cerr << "Problem with NLL location!" << std::endl;
    return false;
  }
  
  updatePicksWithNLL(nllPickList, this);
  
  updateLocationWithNLL(cart, sol, this);
  
  return true;
}

//--------------------------------------------------

void Event::write(std::ostream& os,
                  const Cartesian& cart,
                  float (*ttFunc)(const Position&, const Position&, float&),
                  float pToSVelocityRatio,
                  std::map<std::string, int>& stationMap,
                  std::list<Position> stationPositionList)
{
  std::list<MyPick>::iterator miter;
  float lat, lon, ratio;
  TIME utl;
  int id, idcount;
  Position stationPosition;
  std::list<Position>::iterator piter;
  bool foundStation = false;
  std::map<std::string, int>::iterator iter;
	float dummy;
	
	std::cerr << "writing an event" << std::endl; // call method to create origin, then pass thru callback

  cart.toLonLat(location.x, location.y, &lon, &lat);
  Position latlon = Position(lat, lon, location.z);
  utl = do2date(originTime);

  os << std::setprecision(2) << std::fixed << originTime << " ";
  
  os << utl.yr << " "
     << std::setw(2) << std::setfill('0') << utl.mo << " "
     << std::setw(2) << std::setfill('0') << utl.day << " "
     << std::setw(2) << std::setfill('0') << utl.hr << " "
     << std::setw(2) << std::setfill('0') << utl.mn << " "
     << std::setw(5) << std::setprecision(2) << std::fixed << std::setfill('0') << utl.sec << " ";
  
  os.unsetf(std::ios_base::fixed);
  os << std::setprecision(6) << latlon << " "
     << pickList.size() << " "
     << gap << " "
     << rms << std::endl;
  
  for (miter = pickList.begin(); miter != pickList.end(); miter++) {
    if ((iter = stationMap.find(miter->stationCode)) == stationMap.end()) {
      std::cerr << "Event::write: unknown station code found in test pick" << std::endl;
      return;
    }
    else {
      id = iter->second;
    }
    foundStation = false;
    for (piter = stationPositionList.begin(), idcount = 0; piter != stationPositionList.end(); piter++, idcount++) {
      if (id == idcount) {
        foundStation = true;
        stationPosition = *piter;
        break;
      }
    }
    if (!foundStation) {
      std::cerr << "Event::write: cannot find the station position" << std::endl;
      return;
    }
    if (miter->isSPhase()) {
      ratio = pToSVelocityRatio;
    }
    else {
      ratio = 1.0f;
    }
    miter->writeAndreas(os, ratio * ttFunc(location, stationPosition, dummy));
    os << std::endl;
  }
}

//--------------------------------------------------

#ifdef MY_DEBUG
void Event::writeDebug(std::ostream& os, const Cartesian& cart)
{
  std::list<MyPick>::iterator miter;
  float lat, lon;
  
  cart.toLonLat(location.x, location.y, &lon, &lat);
  Position latlon = Position(lat, lon, location.z);
 
  os << "Origin: " << originTime << std::endl;
  os << "Location: " << latlon << std::endl;
  os << "Error: " << error << std::endl;
  os << pickList.size() << " Picks: ";
  for (miter = pickList.begin(); miter != pickList.end(); miter++) {
    miter->writeDebug(os);
    os << " ";
  }
  os << std::endl;
}
#endif