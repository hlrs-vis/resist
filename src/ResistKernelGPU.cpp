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
#include <omp.h>

using namespace std;


#ifdef HUNTER
#pragma omp requires unified_shared_memory
#endif




// Structure to represent a cell in the grid
struct Cell {
	int x, y;
	int cost;
};

// Define the node structure for the linked list
struct Node {
	Cell data;
	Node* next = nullptr;
};

class Stack {
private:
	// This variable keeps track of the stack size
	static const int numNodes = 4000;

public:
	// Top-of-stack pointer
	Node* top = nullptr;
	Node nodes[numNodes];
	int numFree;
	Node* freeNodes[numNodes];
	int size;


	// Constructor to initialize an empty stack
	Stack() {
		top = NULL; // Initialize top pointer to NULL for an empty stack
		size = 0; // Initialize the size of the stack to zero
		for (int i = 0; i < numNodes; i++)
		{
			freeNodes[i] = &nodes[i];
		}
		numFree = numNodes;
	}
	Node* getNode()
	{
		return(freeNodes[--numFree]);
	}
	void returnNode(Node* n)
	{
		freeNodes[numFree++] = n;
	}

	// Function to push an element onto the stack
	void push(Cell x) {
		Node* new_Node = getNode();
		if (new_Node == nullptr)
		{
			//std::cerr << "bad alloc" << std::endl;
			return;
		}
		new_Node->data = x;
		new_Node->next = nullptr;
		Node* insertpoint = top;
		Node* lastinsertpoint = nullptr;
		while (insertpoint && (x.cost > insertpoint->data.cost))
		{
			if (insertpoint->next == nullptr)
			{
				new_Node->next = nullptr;
				insertpoint->next = new_Node;
				size++;
				return;
			}
			lastinsertpoint = insertpoint;
			insertpoint = insertpoint->next;
		}
		if (insertpoint != nullptr)
		{
			new_Node->next = insertpoint;
			if (insertpoint == top)
			{
				top = new_Node;
			}
			if (lastinsertpoint)
			{
				lastinsertpoint->next = new_Node;
			}
			size++;
			return;
		}
		new_Node->next = top;
		top = new_Node;
		size++;
	}


	// Function to pop an element from the stack
	void pop() {
		if (top == NULL) {
			//cout << "Stack is empty!" << endl; // Display message if the stack is empty
			return;
		}
		Node* temp = top; // Store the current top node
		top = top->next; // Move top to the next node
		size--; // Decrement the size of the stack
		returnNode(temp); // Delete the previous top node
	}


};

// Function to check if a cell is within grid bounds
inline bool isValid(int x, int y, int rows, int cols) {
	return x >= 0 && x < cols && y >= 0 && y < rows;
}


// Function to compute the minimum cost distance matrix sum

int getcostsum(int* costGrid, int rows, int cols, int startX, int startY, int cellSi, int sigma)
{
	// Directions for 8-neighbor movement
	const int dx[8] = { -1, -1, -1, 0, 0, 1, 1, 1 };
	const int dy[8] = { -1, 0, 1, -1, 1, -1, 0, 1 };
	const double moveCost[8] = { sqrt(2), 1, sqrt(2), 1, 1, sqrt(2), 1, sqrt(2) };

	// Initialize the distance matrix with infinity
	int* minCost = new int[rows * cols];

#ifdef HUNTER
#pragma omp for collapse(2)
#endif
	for (int co = 0; co < cols; co++)
		for (int ro = 0; ro < rows; ro++)
		{
			minCost[co * rows + ro] = 100000;
		}

	minCost[startX * rows + startY] = costGrid[startX * rows + startY];

	// Min-heap priority queue to store cells to be processed
	Stack pq;
	pq.push({ startX, startY, costGrid[startX * rows + startY] });

	while (pq.top) {
		Cell current = pq.top->data;
		pq.pop();

		int x = current.x;
		int y = current.y;
		double currentCost = current.cost;

		// If the current cost is greater than the recorded minimum cost, skip processing
		if (currentCost > minCost[x * rows + y]) continue;

		// Explore all 8 neighbors
		for (int i = 0; i < 8; ++i)
		{
			int nx = x + dx[i];
			int ny = y + dy[i];

			if (isValid(nx, ny, rows, cols)) {
				// Calculate the movement cost
				int newCost = minCost[x * rows + y] + costGrid[nx * rows + ny] * moveCost[i];

				// If a lower cost path is found, update and push to the priority queue
				if (newCost < minCost[nx * rows + ny])
				{
					minCost[nx * rows + ny] = newCost;
					pq.push({ nx, ny, newCost });
				}
			}
		}
	}

	//Summe der Widerstände cs berechnen
	double cs = 0;
	int t2 = 2*((sigma) * (sigma));
#ifdef HUNTER
#pragma omp for collapse(2)
#endif
	for (int ro = 0; ro < rows; ro++)
		for (int co = 0; co < cols; co++)
		{

			double A_x, A_y, B_x, B_y;
			A_x = startX * cellSi;
			A_y = startY * cellSi;
			B_x = co * cellSi;
			B_y = ro * cellSi;
			double Dx = A_x - B_x;
			double Dy = A_y - B_y;
			double length2 = (Dx * Dx + Dy * Dy);
			double weight = exp(-(length2) / t2);
			if (length2 <= t2)
			{
#ifdef HUNTER
#pragma omp atomic
#endif
				cs = cs + (weight * minCost[co * rows + ro]);
			}

		}
	delete[] minCost;
	return (int)cs;
}


// MAIN--------------------------------------------------------------------------
// MAIN--------------------------------------------------------------------------
// 
int main(int argc, char** argv)
{
	int row, zeilength;
	int p, done;
	string line, anfang, str1, str2;


	//INITIALISIERUNGEN


	if (argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " InputFile.txt" << std::endl;
		return -1;
	}

	ifstream infile;
	stringstream inStr;
	infile.open(argv[1]);

	ofstream outfile;
	outfile.open(std::string(argv[1]) + "resist.asc");



	std::cout.precision(10);
	outfile.precision(10);

	std::cout << "RsistKernel" << endl;
	std::cout << "..\\... .txt -> ..\\.... .txt" << endl;
	std::cout << "Author/(C): svr@ilpoe.uni-stuttgart.de" << endl;
	std::cout << "No guarantee against unexpected runtime exceptions, usage at your own risk." << endl;
	std::cout << "No responsibility for validity of results." << endl << endl;

	int r = 800;
	int m = 2;
	int sigma = r/m;

	std::cout << "Standard-Deviation (m): " << sigma << endl;
	std::cout << "Multiplier: " << m << endl;
	std::cout << "Kernel-radius: " << r << endl;

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

	int* grid = new int[nrows * ncols];

	//ASCII-GRID LESEN

	std::cout << "Reading Raster..." << endl;
	string subs;
	int  aNum;
	size_t pos, posIs;

	zeilength = ncols * 2 + 1;
	for (row = 0; row < nrows; row++)
	{
		getline(infile, line);
		posIs = 0;
		for (int i = 0; i < ncols; i++)
		{
			pos = line.find(' ', posIs);
			subs = line.substr(posIs, pos);
			aNum = stoi(subs);
			if (aNum == -9999)
			{
				aNum = 9999;
			}
			grid[i * nrows + row] = aNum;
			posIs = pos + 1;
		}
		done = row * 100 / nrows;
		std::cout << "\r" << done << "%";
	}

	std::cout << "done." << endl << "Processing..." << endl;



	//SCHLEIFE ÜBER ALLE ZELLEN 	
	int step = 1;
	int* costsums = new int[ncols];

	p = (r / cellsize) - 1;


	for (int rowA = 0; rowA < nrows; rowA++)
	{
		double t = 0.0;
		double tb, te;

		tb = omp_get_wtime();
#ifdef HUNTER
#pragma omp target teams shared(grid, costsums)      
#pragma omp distribute 
#else
#ifdef HAWK
#pragma omp parallel
#pragma omp for
#endif
#endif
		for (int colA = 0; colA < ncols; colA++)
		{

			//UMGEBUNGSMATRIX HERSTELLEN
			int* umge = new int[2 * p * 2 * p];

			int rmin = max(0, rowA - p);
			int rmax = min(nrows, rowA + p);
			int cmin = max(0, colA - p);
			int cmax = min(ncols, colA + p);
			int sr = rmax - rmin;
			int sc = cmax - cmin;

#ifdef HUNTER
#pragma omp parallel for collapse(2)
#endif
			for (int colB = cmin; colB < cmax; colB++)
				for (int rowB = rmin; rowB < rmax; rowB++)
				{
					umge[(colB - cmin) * sr + (rowB - rmin)] = grid[colB * nrows + rowB];
				}

			//BERECHNE UND SPEICHERE KOSTENSUMME
			costsums[colA] = getcostsum(umge, sr, sc, colA - cmin, rowA - rmin, cellsize, sigma);
			delete[] umge;


		}

		te = omp_get_wtime();
		t = te - tb;
		std::cout << "\n" << "Time of kernel: " << t << endl;

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

	infile.close();
	outfile.close();
	std::cout << "\r" << "100% done.";
}

