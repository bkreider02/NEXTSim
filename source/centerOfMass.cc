#include <iostream>
#include <fstream>

#include "G4Step.hh"

#include "nDetDetector.hh"
#include "centerOfMass.hh"
#include "vertilon.hh"

const double coeff = 1.23984193E-3; // hc = Mev * nm

///////////////////////////////////////////////////////////////////////////////
// class centerOfMass
///////////////////////////////////////////////////////////////////////////////

centerOfMass::~centerOfMass(){
}

centerOfMass centerOfMass::clone() const {
	centerOfMass retval;
	retval.response = response.clone();
	for(size_t i = 0; i < 4; i++)
		retval.anodeResponse[i] = anodeResponse[i].clone();
	if (getNumColumns() == 8 && getNumRows() == 8) {
		for (size_t i = 0; i < 8; i++) {
			for (size_t j = 0; j < 8; j++) {
				retval.pixelResponse[i][j] = pixelResponse[i][j].clone();
			}
		}
	}
	retval.gainMatrix = gainMatrix;
	retval.countMatrix = countMatrix;
	return retval;
}

G4ThreeVector centerOfMass::getCenter() const {
	return (totalMass > 0 ? (1/totalMass)*center : center);
}

double centerOfMass::getCenterX() const {
	return (totalMass > 0 ? (1/totalMass)*center.getX() : 0);
}

double centerOfMass::getCenterY() const {
	return (totalMass > 0 ? (1/totalMass)*center.getY() : 0);
}

double centerOfMass::getCenterZ() const {
	return (totalMass > 0 ? (1/totalMass)*center.getZ() : 0);
}

bool centerOfMass::getCenterSegment(G4ThreeVector &pos, short &col, short &row) const {

	if (pmtHasGaps) {
		double lowerBoundX, upperBoundX;
		double lowerBoundY, upperBoundY;

		double xpos = pos.getX() + activeWidth/2;
		double ypos = pos.getY() + activeHeight/2;

		lowerBoundX = 0;
		upperBoundX = pixelWidth;

		// column and row set to -1 (not touching the PMT) by default
		col = -1;
		row = -1;

		// determine (if any) column/row
		for (int i=0; i<Ncol; i++) {
			if ((xpos > lowerBoundX) && (xpos < upperBoundX)) {
				col = i;
			}

			// reset Y bounds
			lowerBoundY = 0;
			upperBoundY = pixelHeight;
			for (int j=0; j<Nrow; j++) {
				if ((ypos > lowerBoundY) && (ypos < upperBoundY)) {
					row = j;
				}
				lowerBoundY = lowerBoundY + pmtGapThickness + pixelHeight;
				upperBoundY = upperBoundY + pmtGapThickness + pixelHeight;
			}
			if (col != -1)
				break;
			
			lowerBoundX = lowerBoundX + pmtGapThickness + pixelWidth;
			upperBoundX = upperBoundX + pmtGapThickness + pixelWidth;
		}

		return (col >= 0 && col < Ncol && row >= 0 && row < Nrow);
	}
	else {	
		double xpos = (pos.getX()+activeWidth/2)/pixelWidth;
		double ypos = (pos.getY()+activeHeight/2)/pixelHeight;
		
		col = (short)floor(xpos);
		row = (short)floor(ypos);
		
		//std::cout << " activeWidth=" << activeWidth << ", activeHeight=" << activeHeight << std::endl;
		//std::cout << " pixelWidth=" << pixelWidth << ", pixelHeight=" << pixelHeight << std::endl;
		//std::cout << " x=" << xpos << ", y=" << ypos << ", col=" << col << ", row=" << row << std::endl;
		
		return ((col >= 0 && col < Ncol) && (row >= 0 && row < Nrow));
	}
}

bool centerOfMass::getCenterSegment(short &col, short &row) const {
	G4ThreeVector pos = this->getCenter();
	return this->getCenterSegment(pos, col, row);
}	

void centerOfMass::getAnodeCurrents(double *array) const {
	for(size_t i = 0; i < 4; i++){
		array[i] = anodeCurrent[i];
	}
}

double centerOfMass::getReconstructedX() const {
	return ((anodeCurrent[0]+anodeCurrent[1])-(anodeCurrent[2]+anodeCurrent[3]))/(anodeCurrent[0]+anodeCurrent[1]+anodeCurrent[2]+anodeCurrent[3]);
}

double centerOfMass::getReconstructedY() const {
	return ((anodeCurrent[1]+anodeCurrent[2])-(anodeCurrent[3]+anodeCurrent[0]))/(anodeCurrent[0]+anodeCurrent[1]+anodeCurrent[2]+anodeCurrent[3]);
}

double centerOfMass::getSipmX() const {
	double total = 0,x=0;
	for(int i=0;i<Ncol;i++){
		for(int j=0;j<Nrow;j++){
			total += (Ncol-1)*countMatrix[i][j];
			x += i*countMatrix[i][j];
		}
	}
	return x/total;
}

double centerOfMass::getSipmY() const {
	double total = 0,y=0;
	for(int i=0;i<Ncol;i++){
		for(int j=0;j<Nrow;j++){
			total += (Nrow-1)*countMatrix[i][j];
			y += (Nrow-1-j)*countMatrix[i][j];
		}
	}
	return y/total;
}

short centerOfMass::setNumColumns(const short &col_){ 
	Ncol = col_;
	if (pmtHasGaps)
		pixelWidth = (activeWidth - pmtGapThickness*(Ncol-1))/Ncol;
	else
		pixelWidth = activeWidth / (Ncol > 0 ? Ncol : 1);
	return Ncol; 
}

short centerOfMass::setNumRows(const short &row_){
	Nrow = row_;
	if (pmtHasGaps)
		pixelHeight = (activeHeight - pmtGapThickness*(Nrow-1))/Nrow;
	else
		pixelHeight = activeHeight / (Nrow > 0 ? Nrow : 1);
	return Nrow; 
}

double centerOfMass::setActiveAreaWidth(const double &width_){
	activeWidth = width_;
	if (pmtHasGaps) 
		pixelWidth = (activeWidth - pmtGapThickness*(Ncol-1))/Ncol;
	else
		pixelWidth = activeWidth / (Ncol > 0 ? Ncol : 1);
	return activeWidth;
}

double centerOfMass::setActiveAreaHeight(const double &height_){
	activeHeight = height_;
	if (pmtHasGaps) 
		pixelHeight = (activeHeight - pmtGapThickness*(Nrow-1))/Nrow;
	else
		pixelHeight = activeHeight / (Nrow > 0 ? Nrow : 1);
	return activeHeight;
}

void centerOfMass::setSegmentedPmt(const nDetDetectorParams *params){
	if(getNumColumns() <= 0 || getNumRows() <= 0)
		return;

	Ncol = params->GetNumPmtColumns();
	Nrow = params->GetNumPmtRows();
	
	if (params->PmtHasGaps()) {
		pmtHasGaps = true;
		activeWidth = params->GetPmtWidth();
		activeHeight = params->GetPmtHeight();
		pixelWidth = params->GetPmtPixelWidth();
		pixelHeight = params->GetPmtPixelHeight();
		pmtGapThickness = params->GetPmtGapThickness();
	}
	else {
		activeWidth = params->GetPmtWidth();
		activeHeight = params->GetPmtHeight();
		pixelWidth = activeWidth / Ncol;
		pixelHeight = activeHeight / Nrow;
	}
	
	// Setup the anode gain matrix.
	gainMatrix.clear();
	// std::cout<<"Clearing count matrix"<<std::endl;
	// countMatrix.clear();
	for(short i = 0; i < Ncol; i++){
		gainMatrix.push_back(std::vector<double>(Nrow, 100));
		countMatrix.push_back(std::vector<int>(Nrow, 0));
	}
}

bool centerOfMass::loadSpectralResponse(const char *fname){
	return response.loadSpectralResponse(fname);
}

bool centerOfMass::loadGainMatrix(const char *fname){
	if(gainMatrix.empty() || Ncol*Nrow == 0) {/*std::cout<<"Gain Matrix Empty"<<std::endl;*/ return false; }
	std::ifstream gainFile;
	gainFile.open(fname);
	if(!gainFile.good()) { std::cout<<"Gain File Bad"<<std::endl; return false;}
	double readval;
	for(short col = 0; col < Ncol; col++){
		for(short row = 0; row < Nrow; row++){
			gainFile >> readval;
			if(gainFile.eof()){
				gainFile.close();
				return false;
			}
			gainMatrix[col][row] = readval;
		}
	}
	
	gainFile.close();
	
	return true;
}

void centerOfMass::copySpectralResponse(centerOfMass *other){
	response.copySpectralResponse(other->getPmtResponse()->getSpectralResponse());
}

void centerOfMass::copySpectralResponse(const centerOfMass *other){
	response.copySpectralResponse(other->getConstPmtResponse()->getConstSpectralResponse());
}

void centerOfMass::copyGainMatrix(centerOfMass *other){
	other->getGainMatrix(gainMatrix);
}

void centerOfMass::copyGainMatrix(const centerOfMass *other){
	other->getGainMatrix(gainMatrix);
}

void centerOfMass::clear(){
	Npts = 0;
	NnotDetected = 0;
	tSum = 0;
	lambdaSum = 0;
	totalMass = 0;
	center = G4ThreeVector();
	t0 = std::numeric_limits<double>::max();	
	response.clear();
	for(size_t i = 0; i < 4; i++){
		anodeCurrent[i] = 0;
		anodeResponse[i].clear();
	}
	countMatrix.clear();
	for(short i = 0; i < Ncol; i++){
		countMatrix.push_back(std::vector<int>(Nrow, 0));
	}
}

bool centerOfMass::addPoint(const double &energy, const double &time, const G4ThreeVector &position, const double &mass/*=1*/){
	double wavelength = coeff/energy; // in nm
	if(Ncol < 0 && Nrow < 0){ // Default behavior
		center += mass*position;	
		
		// Add the PMT response to the "digitized" trace
		response.addPhoton(time, wavelength);
		
		// Add the "mass" to the others
		totalMass += mass;		
	}
	else{ // Segmented PMT behavior
		short xpos, ypos;
		G4ThreeVector pos = position;
		//std::cout << "pre: " << pos.getX() << "\t" << pos.getY() << std::endl;
		if(this->getCenterSegment(pos, xpos, ypos)){
			// Get the gain of this anode.
			double gain = getGain(xpos, ypos);
			increment(xpos, ypos);

			// Add the anger logic currents to the anode outputs.
			double *current = getCurrent(xpos, ypos);
			if(current){
				for(size_t i = 0; i < 4; i++){
					anodeCurrent[i] += gain*mass*current[i];
				}
				for(size_t i = 0; i < 4; i++){
					anodeResponse[i].addPhoton(time, wavelength, gain*mass*(current[i]));
				}

				if (getNumRows() == 8 && getNumColumns() == 8) { // only record pixel response for 8x8 pmts

					for (size_t i = 0; i < 8; i++) {
						for (size_t j = 0; j < 8; j++) {
							if (i == xpos && j == ypos) {
								pixelResponse[i][j].addPhoton(time, wavelength, gain*mass);
							}
							//else
							//	pixelResponse[i][j].addPhoton(time, 0, 0); // ignore pixels that aren't hit
						}
					}
				}
			}
			
			// Compute resistor network leakage current. This is unnecessary for symmetric leakage... CRT
			/*const double leakage[3][3] = {{1E-3, 1E-2, 1E-3},
			                              {1E-2, 1.00, 1E-2},
			                              {1E-3, 1E-2, 1E-3}};
			for(short anodeX = -1; anodeX <= 1; anodeX++){
				for(short anodeY = -1; anodeY <= 1; anodeY++){
					double *current = getCurrent(xpos+anodeX, ypos+anodeY);
					if(current){ // Add the anger logic currents to the anode outputs.
						for(size_t i = 0; i < 4; i++){
							anodeCurrent[i] += gain*leakage[anodeX+1][anodeY+1]*current[3-i];
						}
					}
				}
			}*/
			
			// Add the PMT response to the "digitized" trace
			response.addPhoton(time, wavelength, gain);

			// Add the "mass" to the others weighted by the individual anode gain
			center += pos; //add pos calculated based on anode current
			totalMass += mass;
		}
	}
	
	tSum += time;
	lambdaSum += wavelength;
	if(time < t0) t0 = time;

	Npts++;
	
	return true;
}

void centerOfMass::printCounts() const {
	for(short i = Nrow-1; i >= 0; i--){
		for(short j = 0; j < Ncol; j++){
			std::cout << countMatrix[j][i] << "\t";
		}
		std::cout << std::endl;
	}		
}

void centerOfMass::print() const {
	if(!empty()){
		std::cout << "M=" << totalMass << ", c=(" << getCenterX() << ", " << getCenterY() << ", " << getCenterZ() << ")\n";
		std::cout << " t0=" << t0 << ", tAvg=" << tSum/Npts << std::endl;
	}
}

void centerOfMass::increment(const int &x, const int &y){
	if((x < 0 || x >= Ncol) || (y < 0 || y >= Nrow) || countMatrix.empty()) return;
	countMatrix[x][y]++;
}

double centerOfMass::getGain(const int &x, const int &y){
	if((x < 0 || x >= Ncol) || (y < 0 || y >= Nrow) || gainMatrix.empty()) return 0;
	return gainMatrix[x][y]/100;
}

double *centerOfMass::getCurrent(const int &x, const int &y){
	if((x < 0 || x >= 8) || (y < 0 || y >= 8)) return NULL;
	return vertilon::currents[x][y];
}
