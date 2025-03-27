// ResistKernelC.cpp : Diese Datei enthält die Funktion "main". Hier beginnt und endet die Ausführung des Programms.
//

#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h> 
#include <algorithm>
#include <string>
#include <time.h>
#include <sstream>
#include <cmath>

using namespace std;

//#pragma omp requires unified_shared_memory

// HELPERS--------------------------------------------------------------------------
// HELPERS

	

// Structure to represent a cell in the grid
struct Cell {
	int x, y;
	double cost;
	bool operator>(const Cell& other) const {
		return cost > other.cost;
	}
};


// Function to check if a cell is within grid bounds
inline bool isValid(int x, int y, int rows, int cols) {
	return x >= 0 && x < rows && y >= 0 && y < cols;
}


// MUSIK--------------------------------------------------------------------------
// MUSIK
// 
// Function to compute the minimum cost distance matrix sum

int getcostsum(vector<vector<double>>& costGrid, int startX, int startY, int cellSi, int rad, int maxWid)
{
// Directions for 8-neighbor movement
const int dx[8] = { -1, -1, -1, 0, 0, 1, 1, 1 };
const int dy[8] = { -1, 0, 1, -1, 1, -1, 0, 1 };
const double moveCost[8] = { sqrt(2), 1, sqrt(2), 1, 1, sqrt(2), 1, sqrt(2) };
	int rows = costGrid.size();
	int cols = costGrid[0].size();

	// Initialize the distance matrix with infinity
	vector<vector<double>> minCost(rows, vector<double>(cols, numeric_limits<double>::infinity()));
	minCost[startX][startY] = costGrid[startX][startY];

	// Min-heap priority queue to store cells to be processed
	priority_queue<Cell, vector<Cell>, greater<Cell>> pq;
	pq.push({ startX, startY, costGrid[startX][startY] });

	while (!pq.empty()) {
		Cell current = pq.top();
		pq.pop();

		int x = current.x;
		int y = current.y;
		double currentCost = current.cost;

		// If the current cost is greater than the recorded minimum cost, skip processing
		if (currentCost > minCost[x][y]) continue;

		// Explore all 8 neighbors
		for (int i = 0; i < 8; ++i) 
		{
			int nx = x + dx[i];
			int ny = y + dy[i];

			if (isValid(nx, ny, rows, cols)) {
				// Calculate the movement cost
				double newCost = minCost[x][y] + costGrid[nx][ny] * moveCost[i];

				// If a lower cost path is found, update and push to the priority queue
				if (newCost < minCost[nx][ny]) 
				{
					minCost[nx][ny] = newCost;
					pq.push({ nx, ny, newCost });
				}
			}
		}
	}

	//Summe der Widerstände cs berechnen
	int cs = 0;
	for (int ro = 0; ro < rows; ro++)
	{
		for (int co = 0; co < cols; co++)
		{
			double A_x, A_y, B_x, B_y;
			A_x = startX * cellSi;
			A_y = startY * cellSi;
			B_x = co * cellSi;
			B_y = ro * cellSi;
			double Dx = A_x - B_x;
			double Dy = A_y - B_y;
			int d = sqrt(Dx*Dx + Dy*Dy);
			if (d <= rad)
			{
				if (cs <= maxWid)
				{
					cs = cs + minCost[ro][co];
				}
			}
		}
	}	
	return cs;
}


// MAIN--------------------------------------------------------------------------
// MAIN--------------------------------------------------------------------------
// 
int main()
{
	int cs, fd, anzges, anz1, row, i, s, zeilength;
	int p, q, done, costsum;
	double d, dsum, sumges, A_x, A_y, B_x, B_y, xout, yout, wcc, dis;
	string line, anfang, str1, str2;
	
	
	//INITIALISIERUNGEN

	//int(*grid)[32767] = new int[32767][32767];
	//int(*umge)[32767] = new int[32767][32767];

	//const vector<vector<double>>& grid;
	//const vector<vector<double>>& umge;
	


	ifstream infile;
	stringstream inStr;
	infile.open("/zhome/academic/rus/rus/woess/ilpoe/BeispielMatrix.txt");

	ofstream outfile;
	outfile.open("/zhome/academic/rus/rus/woess/ilpoe/WidObfDet30_resist.asc");

	std::cout.precision(10);
	outfile.precision(10);

	std::cout << "RsistKernel" << endl;
	std::cout << "..\\... .txt -> ..\\.... .txt" << endl;
	std::cout << "Author/(C): svr@ilpoe.uni-stuttgart.de" << endl;
	std::cout << "No guarantee against unexpected runtime exceptions, usage at your own risk." << endl;
	std::cout << "No responsibility for validity of results." << endl << endl;

	int r;
	int reslim;
	r = 1000;
	reslim = 1000;
	std::cout << "Kernel-radius: " << r << endl;
	std::cout << "Max-Resist: " << reslim << endl;

	//ASCII-HEADER LESEN	
	getline(infile, line);
	anfang = line.substr(0, 13);
	line.erase(0, 13);
	int ncols = atoi(line.c_str());
	std::cout << "ncols:" << ncols << endl;
	outfile << anfang << " " << ncols << endl;

	getline(infile, line);
	anfang = line.substr(0, 13);
	line.erase(0, 13);
	int nrows = atoi(line.c_str());
	std::cout << "nrows:" << nrows << endl;
	outfile << anfang << " " << nrows << endl;

	getline(infile, line);
	anfang = line.substr(0, 13);
	line.erase(0, 13);
	double xll = atoi(line.c_str());
	std::cout << "xll:" << xll << endl;
	outfile << anfang << " " << xll << endl;

	getline(infile, line);
	anfang = line.substr(0, 13);
	line.erase(0, 13);
	double yll = atoi(line.c_str());
	std::cout << "yll: " << yll << endl;
	outfile << anfang << " " << yll << endl;

	getline(infile, line);
	anfang = line.substr(0, 13);
	line.erase(0, 13);
	int cellsize = atoi(line.c_str());
	std::cout << "cellsize: " << cellsize << endl;
	outfile << anfang << " " << cellsize << endl;

	getline(infile, line);
	outfile << line << endl;

	vector<vector<int>> grid(nrows, vector<int>(ncols, numeric_limits<int>::infinity()));
	
	//ASCII-GRID LESEN

	std::cout << "Reading Raster..." << endl;
	string subs;
	int pos, posIs, aNum;

	zeilength = ncols * 2 + 1;
	for (row = 0; row < nrows; row++)
	{
		getline(infile, line);			
		for (i = 0; i < ncols; i++)
		{
			posIs = 0;
			pos = line.find(' ');
			subs = line.substr(posIs, pos);				
			aNum = stoi(subs);
			if (aNum == -9999)
			{
				aNum = 9999;
			}
			grid[row][i] = aNum;
			posIs = pos;
		}
		done = row * 100 / nrows;
		std::cout << "\r" << done << "%";
	}

	std::cout << "done." << endl << "Processing..." << endl;
	
	

	//SCHLEIFE ÜBER ALLE ZELLEN 	
	int step = 1;
	vector<int> costsums(ncols, 0);

	p = (r / cellsize) - 1;


	for (int rowA = 0; rowA < nrows; rowA++)
	{
#pragma omp parallel
#pragma omp target teams distribute
		for (int colA = 0; colA < ncols; colA++)
		{

			//UMGEBUNGSMATRIX HERSTELLEN
			vector<vector<double>> umge(2 * p, vector<double>(2 * p, numeric_limits<int>::infinity()));
			int v = 0;

			int rmin = max(0, rowA - p);
			int rmax = min(nrows, rowA + p);
			int cmin = max(0, colA - p);
			int cmax = min(ncols, colA + p);

			for (int rowB = rmin; rowB < rmax; rowB++)
			{
				int w = 0;
				for (int colB = cmin; colB < cmax; colB++)
				{					
					umge[v][w] = grid[rowB][colB];
					w++;
				}
				v++;
			}

			//BERECHNE UND SPEICHERE KOSTENSUMME

			costsums[colA] = getcostsum(umge, colA - cmin, rowA - rmin, cellsize, r, reslim);
		
		
		}

		for (int colA = 0; colA < ncols; colA++)
		{
			//done = step * 100 / (nrows*ncols);
			step++;

			outfile << costsums[colA];
			if (colA < ncols)
			{
				outfile << " ";
			}

		}
		std::cout << "\n" << step << " (" << nrows * ncols << ")";
		std::cout << "\n" << rowA << "/" << nrows;
		outfile << endl;

		
	}

	std::cout << "\r" << "100% done.";
	infile.close();
	outfile.close();
	std::cout << endl << "Enter q to quit:";
	std::cin >> q;
}

