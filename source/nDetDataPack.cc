#include <iostream>

#include "TBranch.h"

#include "nDetDataPack.hh"

void nDetDataPack::setDataAddresses(nDetEventStructure *evt, nDetOutputStructure *out, nDetImplantOutputStructure *outImp, nDetMultiOutputStructure *mult, nDetDebugStructure *debug, nDetTraceStructure *trace){
	evtData = evt;
	outData = out;
	outImplantData = outImp;
	multData = mult;
	debugData = debug;
	traceData = trace;
}

void nDetDataPack::copyData(nDetEventStructure *evt, nDetOutputStructure *out, nDetImplantOutputStructure *outImp, nDetMultiOutputStructure *mult, nDetDebugStructure *debug, nDetTraceStructure *trace) const {
	(*evt) = (*evtData);
	(*out) = (*outData);
	(*outImp) = (*outImplantData);
	(*mult) = (*multData);
	(*debug) = (*debugData);
	(*trace) = (*traceData);
}

bool nDetDataPack::goodEvent() const {
	return (evtData->goodEvent || (debugData->nPhotons[0] > 0 || debugData->nPhotons[1] > 0));
}

int nDetDataPack::getEventID() const {
	return evtData->eventID;
}

void nDetDataPack::clear(){
	evtData->Zero();
	outData->Zero();
	outImplantData->Zero();
	multData->Zero();
	debugData->Zero();
	traceData->Zero();
}
