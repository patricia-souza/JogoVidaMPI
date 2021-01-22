#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include "mpi.h"

using std::cout;
using std::cin;
using std::endl;
using namespace MPI;

const int TAG_HALO_SWAP = 1;
const int squareSize = 2;

inline int getIndex(int, int, int);
void initGrid(int, int, int, bool[]);
int getNeighbors(int, int, int, bool[]);
void step(int, int, bool[], bool[]);
void swapHalo(Cartcomm, Datatype, int, int, bool[]);
void copyGrid(int, int, bool[], bool[]);
void makeImage(int, int, bool[]);
void printMatrix(int, bool[]);
int countCells(int, bool[]);

inline int getIndex(int width, int y, int x)
{
	return x + (width + 2) * y;
}

void initGrid(int height, int width, int probability, bool grid[])
{
	int y, x;
	for(y = 1; y <= height; y++)
	{
		for(x = 1; x <= width; x++)
		{
			int random =  (100.0 * rand()) / RAND_MAX;
			if(random < probability)
				grid[getIndex(width, y, x)] = true;
			else
				grid[getIndex(width, y, x)] = false;
		}
	}
}

int getNeighbors(int width, int y, int x, bool grid[])
{
	int count = 0;
	if(grid[getIndex(width, y - 1, x - 1)])
		count++;
	if(grid[getIndex(width, y - 1, x)])
		count++;
	if(grid[getIndex(width, y - 1, x + 1)])
		count++;
	if(grid[getIndex(width, y, x - 1)])
		count++;
	if(grid[getIndex(width, y, x + 1)])
		count++;
	if(grid[getIndex(width, y + 1, x - 1)])
		count++;
	if(grid[getIndex(width, y + 1, x)])
		count++;
	if(grid[getIndex(width, y + 1, x + 1)])
		count++;
	return count;
}

void step(int height, int width, bool grid[], bool nextGrid[])
{
	int x, y;
	for(y = 1; y <= height; y++)
	{
		for(x = 1; x <= width; x++)
		{
			int index = getIndex(width, y, x);
			int neighbors = getNeighbors(width, y, x, grid);
			if(!grid[index])
			{
				if(neighbors == 3)
					nextGrid[index] = true;
			}
			else if(neighbors < 2 || neighbors > 3)
				nextGrid[index] = false;
		}
	}
}

void swapHalo(Cartcomm cartesian, Datatype rowType, int height, int width, bool grid[])
{
	int top, bottom;
	cartesian.Shift(0, 1, top, bottom);

	//envia linha ghost superior
	grid[getIndex(width, 1, 0)] = grid[getIndex(width, height, width)];
	grid[getIndex(width, 1, width + 1)] = grid[getIndex(width, height, 1)];
	cartesian.Sendrecv(&grid[getIndex(width, 1, 0)], 1, rowType, top, TAG_HALO_SWAP, &grid[getIndex(width, height + 1, 0)], 1, rowType, bottom, TAG_HALO_SWAP);

	//envia linha ghost inferior
	grid[getIndex(width, height, 0)] = grid[getIndex(width, 1, width)];
	grid[getIndex(width, height, width + 1)] = grid[getIndex(width, 1, 1)];
	cartesian.Sendrecv(&grid[getIndex(width, height, 0)], 1, rowType, bottom, TAG_HALO_SWAP, grid, 1, rowType, top, TAG_HALO_SWAP);

	//para as colunas ghosts
	for(int y = 1; y <= height; y++)
	{
		grid[getIndex(width, y, 0)] = grid[getIndex(width, y, width)];
		grid[getIndex(width, y, width + 1)] = grid[getIndex(width, y, 1)];
	}
}

void copyGrid(int height, int width, bool source[], bool target[])
{
	for(int i = 0; i < (width + 2) * (height + 2); i++){
		target[i] = source[i];
	}
}

void printMatrix(int gridSize, bool grid[]){
	for(int y = 1; y <= gridSize; y++){
		for(int x = 1; x <= gridSize; x++){
			cout << grid[getIndex(gridSize,y,x)] << " ";
		}
		cout << endl;
	}
}

int countCells(int gridSize, bool grid[]){
	int cells = 0;
	for(int y = 1; y <= gridSize; y++){
		for(int x = 1; x <= gridSize; x++){
			if(grid[getIndex(gridSize, y, x)] == true) {
				cells++;
			}
		}
	}
	return cells;
}

int main(int argc, char* argv[])
{
	Init(argc, argv);

	int cpus = COMM_WORLD.Get_size();
	int rank = COMM_WORLD.Get_rank();
	srand(time(NULL) + rank);
	int dims[1];
	dims[0] = cpus;
	bool periods[] = {true};
	Cartcomm cartesian = COMM_WORLD.Create_cart(1, dims, periods, false);

	int gridSize;
	int probability;
	int steps;
	int cells; 

	if(rank == 0)
	{
		cout << "Entre com o tamanho do tabuleiro: ";
		cin >> gridSize;
		cout << "Probabilidade para preencher o tabuleiro [0 - 100]: ";
		cin >> probability;
		cout << "Número de gerações: ";
		cin >> steps;
	} 

	COMM_WORLD.Bcast(&gridSize,    1, INT, 0);
	COMM_WORLD.Bcast(&probability, 1, INT, 0);
	COMM_WORLD.Bcast(&steps,       1, INT, 0);

	int width  = gridSize;
	int height = rank < (gridSize % cpus) ? gridSize / cpus + 1 : gridSize / cpus;

	bool* wholeGrid = new bool[(gridSize + 2) * (gridSize + 2)];
	bool* grid      = new bool[(width + 2) * (height + 2)];
	bool* nextGrid  = new bool[(width + 2) * (height + 2)];
	Datatype rowType = BOOL.Create_contiguous((width + 2));
	rowType.Commit();
	for(int i = 0; i < gridSize * gridSize; i++)
	{
		if(i < (width + 2) * (height + 2))
		{
			grid[i] = false;
			nextGrid[i] = false;
		}
		wholeGrid[i] = false;
	}
	initGrid(height, width, probability, grid); 
	swapHalo(cartesian, rowType, height, width, grid); 

	for(int i = 1; i <= steps; i++)
	{
		step(height, width, grid, nextGrid);
		swapHalo(cartesian, rowType, height, width, nextGrid);
		copyGrid(height, width, nextGrid, grid);
		cartesian.Gather(grid + width + 2, height, rowType, wholeGrid + width + 2, height, rowType, 0);
	}

	if(rank == 0)
	{
		int cells = 0; 
		//printMatrix(gridSize, wholeGrid);
		cells = countCells(gridSize, wholeGrid);
		cout << "CÉLULAS VIVAS: " << cells << endl; 
	}

	rowType.Free();
	delete nextGrid;
	delete grid;
	delete wholeGrid;
	Finalize();
	return 0;
}