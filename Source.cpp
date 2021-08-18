//Optimisation program for distribution hub in UK
/*
Author: James Lockwood
Date: 16-12-2018
Program: 
	Given a CSV file of cities organised in columns, such that each row has the 
	following structure: {Name, Type, Population, Lat, Lon},
	this program will find the ideal location to position a distribution hub.
	The user has the opportunity to select how many hubs they want and whether they would 
	like to consider the ports (supply points into consideration) and their future employess' salary.

	For reducing operating costs of each hub, an optimised hillclimbing algorithm is used. 
	For k hubs, weighted k-means clustering is used, with hillclimbing used to position each hub
	efficiently when the cities have been split into clusters.

	The ports are given in a separate file, GBports.csv, and are simply considered as extra cities
	by the program. Each port has a fraction of the total cities population assigned to it based on 
	its yearly tonnage (see Wiki page of List of busiest ports in Europe). This is down to the 
	assumption that total supplies is related to the total amount of (potential) buyers/customers.
	Depending on the fraction of supplies that your company gets from abroad, the weight of the ports 
	can be modulated in the initial input parameters.

	Salaries can also be factored into the costs of running a hub and are give by a CSV file, GBsalaries.csv,
	where the data comes from https://www.theguardian.com/news/datablog/2011/nov/24/wages-britain-ashe-mapped,
	and was correlated with the known positions of the cities. When iterating, the costs of a hub are multiplied
	by the average yearly salary of the nearest town or city.
Notes:
	Developed in Visual Studio 2017
*/

#include<iostream>	//cout,cin
#include<fstream>	//opening files
#include<vector>	//vector
#include<string>	//string
#include<cmath>		//atan2,sin,cos
#include<chrono>	//high_resolution_clock
#include<iomanip>	//setprecision
#include<cfloat>	//contains certain parameters such as DBL_MAX

#define PI 3.141592654
#define R 6370.99056

using namespace std;
using namespace std::chrono;

//track calls of weighted distance and haversine functions
int iter = 0;
int haversineCount = 0;

//user defined specifications
int k = 5;						//number of hubs
double minstep = pow(2, -15);	//minimum step for the hillclimb algorithm
int attempts = 50;				//number of attempts to find tallest hill
int prec = 1000;				//resolution of random numbers
string filename = "GBplaces.csv"; 
string portname = "GBports.csv";
string salaryname = "GBsalaries.csv";
bool salaryChoice = false;

//degree to radian
double dtor(double degree) {
	return PI / 180 * degree;
}

//Radian to degree
double rtod(double radian) {
	return 180 / PI * radian;
}

//max element of a 1-D vector
double max_element(vector<double> vect) {
	double max = 0;
	for (size_t i = 0; i < vect.size(); i++) {
		if (vect[i] > max) max = vect[i];
	}
	return max;
}

//min element of a 1-D vector
double min_element(vector<double> vect) {
	double min = 0;
	for (size_t i = 0; i < vect.size(); i++) {
		if (vect[i] < min) min = vect[i];
	}
	return min;
}

//random, bounded number with resolution of n
double random_number(double lower, double upper, int n) {
	double r;
	r = lower + (rand() % (n + 1) * (1. / n) * (upper - lower));
	return r;
}

//find crude centre of mass based on lat, lon and population
vector<double> COM(vector<vector<double>> cluster) {
	/*Returns coordinates of centre of mass of a cluster {{xi,yi,pi},...}*/

	double xtot = 0;
	double ytot = 0;
	double ptot = 0;

	for (size_t i = 0; i < cluster.size(); i++) {
		xtot += cluster[i][0] * cluster[i][2];
		ytot += cluster[i][1] * cluster[i][2];
		ptot += cluster[i][2];
	}
	
	return { xtot/ptot,ytot/ptot };
}

//Haversine distance formula
double haversine(double lat1, double lon1, double lat2, double lon2) {
	lat1 = dtor(lat1);
	lat2 = dtor(lat2);
	lon1 = dtor(lon1);
	lon2 = dtor(lon2);

	haversineCount++; //track number of times function is called
	return 2*R* asin(sqrt(pow(sin((lat2 - lat1) / 2), 2) + cos(lat1)*cos(lat2)*pow(sin((lon2 - lon1) / 2), 2)));
	
}

//total weighted distance (factor: population) to all towns and cities
double weightedDist(double lon, double lat, vector<double> xs, vector<double> ys, vector<double> ms) {
	double wD = 0;
	double M = 0;
	//add up all weighted distances and keep track of total population
	for (size_t i = 0; i < xs.size(); i++) {
		wD += haversine(lat, lon, ys[i], xs[i])*ms[i];
		M += ms[i];
	}
	//track number of times function is called
	iter++;
	
	return wD / M ;
}
//nearest city's salary
double nearestSalary(double x, double y, vector<double> xsal, vector<double> ysal, vector<double> psal,bool salaries=salaryChoice) {
	double salary;
	if (salaries) {
		double minDist = DBL_MAX;
		double dist;
		int cidx;
		for (size_t i = 0; i < xsal.size(); i++) {
			dist = haversine(y, x, ysal[i], xsal[i]);
			if (dist < minDist) {
				minDist = dist;
				cidx = i;
			}
		}
		salary = psal[cidx];
	}
	else {
		salary = 1;
	}
	return salary;
}
//nearest city
string nearestCity(double x, double y, vector<double> xs, vector<double> ys,vector<string>cs) {
	double minDist = DBL_MAX;
	double dist;
	int cidx{}; //city index
	for (size_t i = 0; i < xs.size(); i++) {
		dist = haversine(y, x, ys[i], xs[i]);
		if (dist < minDist) {
			minDist = dist;
			cidx = i;
		}
	}
	return cs[cidx];
}
//cost of operations: unaveraged weighted distance
double cost(double lon, double lat, vector<double> xs, vector<double> ys, vector<double> ms) {
	double wD = 0;
	for (size_t i = 0; i < xs.size(); i++) {
		//sum all haversine distances weighted with population
		wD += haversine(lat, lon, ys[i], xs[i])*ms[i];
	}
	//divide by 10^6 to improve readibility of costs
	return wD/pow(10,6);

}

vector<double> hillclimb(vector<double> xs, vector<double> ys, vector<double> ps, vector<double> start,vector<double> xsal, vector<double>ysal, vector<double>psal,double step = 1) {
	/*Hillclimb function to find max/min by incrementing position with increasingly small steps*/

	//coordinates
	vector<double> hilltop;
	
	//initial parameters
	int dx, dy;
	double xpos = start[0];
	double ypos = start[1];

	//memory: avoids remeasuring already evaluated points
	double oldi, oldj;
	double oldvalue, value, newvalue;
	bool changed = false;
	
	//find average cost for start position
	value = weightedDist(xpos, ypos, xs, ys, ps)*nearestSalary(xpos, ypos, xsal, ysal, psal);
	//find lowest average cost
	do {
		//check 8 points surrounnding initial guess for first time using new step value
		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				if (i == 0 && j == 0) {}
				else {
					newvalue = weightedDist(xpos + i * step, ypos + j * step, xs, ys, ps)
						*nearestSalary(xpos+i*step,ypos+j*step,xsal,ysal,psal);
					//update value if new value is smaller
					if (newvalue < value) {
						dx = i;
						dy = j;
						value = newvalue;
						changed = true;
						//cout << xpos + dx * step << " " << ypos + dy * step << " " << value << endl;
					}
				}
			}
		}
		//update position if any lower position has just been found
		if (changed) {
			changed = false;
			xpos += step * dx;
			ypos += step * dy;
			
			//avoid measuring 9 points if lower values still being found with current step value
			do {
				oldvalue = value;
				oldi = dx;
				oldj = dy;
				//new position is top or bottom side of previous 9 point square 
				if (dx == 0) {
					int j = oldj;
					//only measure 3 points on line above/below current position
					for (int i = -1; i <= 1; i++) {
						newvalue = weightedDist(xpos + step * i, ypos + step * j, xs, ys, ps)
							*nearestSalary(xpos + i * step, ypos + j * step, xsal, ysal, psal);
						if (newvalue < value) {
							dx = i;
							dy = j;
							value = newvalue;
							changed = true;
						}
					}
				}
				//new position is left or right side of previous 9 point square
				else if (dy == 0) {
					int i = oldi;
					//only measure 3 points on line left/right current position
					for (int j = -1; j <= 1; j++) {
						newvalue = weightedDist(xpos + step * i, ypos + step * j, xs, ys, ps)
							*nearestSalary(xpos + i * step, ypos + j * step, xsal, ysal, psal);
						if (newvalue < value) {
							dx = i;
							dy = j;
							value = newvalue;
							changed = true;
						}
					}
				}
				//new position is on corner of previous 9 point square
				else {
					//only measure 5 points on L-shape (other 4 previously measured with higher average cost)
					//3 elements of row
					int i = oldi;
					for (int j = -1; j <= 1; j++) {
						newvalue = weightedDist(xpos + step * i, ypos + step * j, xs, ys, ps)
							*nearestSalary(xpos + i * step, ypos + j * step, xsal, ysal, psal);
						if (newvalue < value) {
							dx = i;
							dy = j;
							value = newvalue;
							changed = true;
						}
					}
					int j = oldj;
					//2 elements of column (one measured in previous for loop)
					for (int i = -1; i <= 1; i++) {
						if (i != oldi) {
							newvalue = weightedDist(xpos + step * i, ypos + step * j, xs, ys, ps)
								*nearestSalary(xpos + i * step, ypos + j * step, xsal, ysal, psal);
							if (newvalue < value) {
								dx = i;
								dy = j;
								value = newvalue;
								changed = true;
								//cout << xpos + dx * step << " " << ypos + dy * step << " " << value << endl;
							}
						}
					}
				}
				//update position if lower position found
				if (changed) {
					changed = false;
					xpos += step * dx;
					ypos += step * dy;
				}
			} while (value < oldvalue);
		}
		//halve the step value
		step /= 2;
		changed = false;
	} while (step >= minstep);
	
	//add coords and average cost to return vector
	hilltop.push_back(xpos);
	hilltop.push_back(ypos);
	hilltop.push_back(value);
	return hilltop;
}

vector<double> highesthill(vector<double> xs, vector<double> ys, vector<double> ps, vector<double> xsal, vector<double>ysal, vector<double>psal, int tries = attempts) {
	/*Perform hillclimbing at various points to find global minimum*/

	//track which hill is biggest
	vector<double> biggesthill;
	vector<double> currenthill;
	double peak;

	//find furthest cities and towns
	double xmax = max_element(xs);
	double xmin = min_element(xs);
	double ymax = max_element(ys);
	double ymin = min_element(ys);

	//first random point, constrained to furthest points
	double xrand = random_number(xmin, xmax, prec);
	double yrand = random_number(ymin, ymax, prec);
	
	//climb hill and track peak value
	biggesthill = hillclimb(xs, ys, ps,{ xrand,yrand }, xsal, ysal, psal);
	peak = biggesthill[2];
	
	//attempt another n-1 hills (first attempt above) 
	for (int i = 0; i < tries-1; i++) {
		xrand = random_number(xmin, xmax, prec);
		yrand = random_number(ymin, ymax, prec);
		currenthill = hillclimb(xs, ys, ps,{ xrand ,yrand}, xsal, ysal, psal);
		if (currenthill[2] < peak) {
			biggesthill = currenthill;
			peak = biggesthill[2];
		}
	}
	//returns coordinates and average cost for biggest hill
	return biggesthill;
}

vector<vector<double>> kmeans(int k, vector<double> xs, vector<double> ys, vector<double> ps, vector<double> xsal, vector<double>ysal, vector<double>psal, int tries = attempts) {
	/*Function that returns vector of coordinates of k hubs and a vector of the overall costs for each hub*/

	//vectors for hub coordinates and for clusters
	vector<vector<double>> hubs;
	vector<vector<vector<double>>> finalclusters;

	//copy vectors to remove initial hub possibilities to avoid duplicates
	vector<double> xc = xs;
	vector<double> yc = ys;

	//choose k random start points from cities and towns
	int n= xs.size();
	for (int i = 0; i < k; i++) {
		//random indices
		int xind = rand() % n;
		int yind = rand() % n;
		//add city or town to initial hub list
		hubs.push_back({ xc[xind],yc[yind] });
		//remove city or town from possible options to ensure no du
		xc.erase(xc.begin() + xind);
		yc.erase(yc.begin() + yind);
		n--;
	}

	//check if cluster has changed
	bool changes = false;
	//track previous cluster populations
	vector<int> memory(k,0);

	do {
		changes = false;
		double dist;
		int cidx;
		//(re)initialise cluster vector and populate with empty vectors
		vector<vector<vector<double>>> clusters;
		for (int cluster = 0; cluster < k; cluster++) {
			clusters.push_back({});
		}

		//assign each node to cluster with nearest centroid
		for (size_t e = 0; e < xs.size(); e++) {
			double mindist = DBL_MAX;
			for (int c = 0; c < k; c++) {
				dist = haversine(ys[e], xs[e], hubs[c][1], hubs[c][0])*ps[e];
				//check if current hub is closer
				if (dist < mindist) {
					mindist = dist;
					cidx = c;
				}
			}
			//add city or town to cluster of nearest hub
			clusters[cidx].push_back({ xs[e],ys[e],ps[e] });

		}
		//update clusters
		finalclusters = clusters;
		//reevaluate hub
		for (int c = 0; c < k; c++) {
			//evaluate centre of mass of each cluster and assign hub to it
			hubs[c] = COM(clusters[c]);
			//track if change in cluster numbers
			if (memory[c] != clusters[c].size()) {
				memory[c] = clusters[c].size();
				changes = true;
			}
		}
	} while (changes);

	//display cluster information
	cout << "The " << k << " clusters have the following populations: ";
	for (size_t i = 0; i < finalclusters.size(); i++) {
		cout<<finalclusters[i].size() <<" ";
	}
	cout << endl;


	//convert clusters from vector of points to {xs,ys,ps} for hillclimbing optimisation
	vector<vector<vector<double>>> hubcoords;
	for (int c = 0; c < k; c++) {
		vector<double> xtemp, ytemp, ptemp;
		for (size_t i = 0; i < finalclusters[c].size(); i++) {
			xtemp.push_back(finalclusters[c][i][0]);
			ytemp.push_back(finalclusters[c][i][1]);
			ptemp.push_back(finalclusters[c][i][2]);
		}
		hubcoords.push_back({ xtemp,ytemp,ptemp });
	}

	//optimise hub placement with hillclimb and track overall costs
	vector<double> costs;
	for (int c = 0; c < k; c++) {
		hubs[c] = hillclimb(hubcoords[c][0], hubcoords[c][1], hubcoords[c][2], { hubs[c][0],hubs[c][1] },xsal,ysal,psal);
		costs.push_back(cost(hubs[c][0], hubs[c][1], hubcoords[c][0], hubcoords[c][1], hubcoords[c][2]));
	}
	hubs.push_back(costs);

	return hubs;
}

int ReadData(string filename, vector<double> & xs, vector<double> & ys, vector<double> & ps,vector<string> & cs,double weight=1) {
	ifstream file(filename);

	string line, ptemp, ct;

	double xt, yt, pt;

	if (file.is_open()) {
		getline(file, line);
		while (!file.eof()) {
			getline(file, line);
			int end = line.find(',', 0);
			vector<string> elements;
			for (int i = 0; i < 5; i++) {
				elements.push_back(line.substr(0, end));
				line = line.substr(end + 1, line.length() - end - 1);
				end = line.find(',', 0);

			}
			ct = elements[0];
			xt = atof(elements[4].c_str());
			yt = stod(elements[3].c_str());
			pt = stod(elements[2].c_str())*weight;

			//cout << typeid(xt).name() << endl;
			xs.push_back(xt);
			ys.push_back(yt);
			ps.push_back(pt);
			cs.push_back(ct);
		}
		file.close();
		return 0;
	}
	else {
		cout << filename << " did not open. Exiting program" << endl;
		system("pause");
		return 1;
	}
}

int main() {
	//set general program parameters
	cout << fixed;
	cout << setprecision(8);
	srand(time(NULL));

	//Welcome user and ask for inputs
	string ynans; //yes no answer
	string portans; //choose to consider ports
	string weight;
	string salaries; //choose to consider salaries of employees

	cout << "Welcome to the FindADistributionHub-inator2000!\nFind a hub for "
		<< filename << " with " << k << " hubs: would you like to modify any parameters (y/n)" << endl;
	while (ynans != "y" && ynans != "n") {
		cout << "(please ensure answer is valid) ";
		cin >> ynans;
	}
	if (ynans == "y") {
		cout << "Please note we accept no liability for invalid inputs (or uneconomical hub placements for that matter).\n";
		cout << "Number of hubs: ";
		cin >> k;
		cout << "Take major ports into consideration (y/n): ";
		cin >> portans;
		//if user includes ports, ask for weightage
		if (portans == "y") {
			cout << "What fraction of your supply comes from UK ports (double 0-1): ";
			cin >> weight;
		}
		cout << "Take regional salaries into consideration (y/n): ";
		cin >> salaries;
	}
	//data vectors
	vector<double> xs, ys, ps;
	vector<double> xsal, ysal, psal;
	vector<string> cs;
	vector<string> csal;

	//start the clock
	auto start = high_resolution_clock::now();
	//read in the data to the vectors for cities, ports, and salaries
	int readin = ReadData(filename, xs, ys, ps, cs);
	if (portans == "y") {
		int readin2 = ReadData(portname, xs, ys, ps, cs, stod(weight.c_str()));
		if (readin2) exit(1);
	}
	if (readin) {
		//exits program if error opening file
		exit(1);
	}
	if (salaries == "y") {
		salaryChoice = true;
		int readin3 = ReadData(salaryname, xsal, ysal, psal, csal);
		if (readin3) exit(1);
	}

	//if user selects 1 hub
	if (k == 1) {
		//perform optimised hill climb at multiple random points
		vector<double> hub = highesthill(xs, ys, ps, xsal, ysal, psal);
		//display
		cout << "For the specified parameters, with 1 hub, the coordinates are " <<
			hub[1] << " (lat) and " << hub[0] << " (lon) with an average cost of " << hub[2] << endl;
		cout << "Nearest city to hub: " << nearestCity(hub[0], hub[1], xs, ys, cs) << endl;
		cout << "The operating costs amount to: " << cost(hub[0], hub[1], xs, ys, ps) << endl;
		cout << "\nHave you considered developing multiple hubs? This could reduce your overall operating costs." << endl;
	}
	//find more than one national hub
	else {
		//find k hubs
		vector<vector<double>> hubs = kmeans(k, xs, ys, ps, xsal, ysal, psal);
		//display
		cout << "For the specified parameters, with " << k
			<< " hubs, the hubs were determined as (lat, lon, average cost, nearest city): " << endl;
		for (int i = 0; i < k; i++) {
			cout << hubs[i][1] << " " << hubs[i][0] << " " << hubs[i][2] << " " << nearestCity(hubs[i][0], hubs[i][1], xs, ys, cs) << endl;
		}
		double sum = 0;

		//display costs
		for (size_t i = 0; i < hubs[k].size(); i++) {
			sum += hubs[k][i];
		}
		cout << "Operating costs (all hubs): " << sum << endl;
	}
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout << "\nTime taken for program execution: "
		<< duration.count() << " microseconds\n" << endl;
	
	cout << "Calls of weighted distance function: "<<iter << endl;
	cout << "Calls of haversine function: "<<haversineCount << endl;
	system("pause");
}
