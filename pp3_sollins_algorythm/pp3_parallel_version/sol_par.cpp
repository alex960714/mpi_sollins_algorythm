#include <stdio.h>
#include "mpi.h"
#include <stack>
#include <fstream>
#include <list>
#include <time.h>
using namespace std;

int *matr;
int *mst;
int *comp;
int *calc_en_vert;
int vert_num, comp_num, f_opt, edges_num;
ofstream os;

int ProcNum, ProcRank;
int *sendcounts, *displs;
int *proc_matr;
int *proc_calc_en_vert;
int proc_vert_num, proc_displs;

double st_time, en_time;
MPI_Status status;

void matr_init()
{
	ifstream is;
	is.open("adj_matr.txt");
	is >> vert_num;
	matr = new int[vert_num*vert_num];

	for (int i = 0; i < vert_num; i++)
		for (int j = 0; j < vert_num; j++)
			is >> matr[i*vert_num + j];

	is.close();
}

void matr_init_rand(int v_num, int power)
{
	srand(time(NULL));
	vert_num = v_num;
	matr = new int[vert_num*vert_num];
	for (int i = 0; i < vert_num; i++)
	{
		matr[i*vert_num+i] = 0;
		for (int j = i + 1; j < vert_num; j++)
		{
			matr[i*vert_num + j] = matr[j*vert_num + i] = rand() % power;
		}
	}

	os << "Adjacency matrix:" << endl;
	for (int i = 0; i < vert_num; i++)
	{
		for (int j = 0; j < vert_num; j++)
			os << matr[i*vert_num + j] << " ";
		os << endl;
	}
	os << "--------------------" << endl;
}

void mem_init()
{
	if (!ProcRank)
	{
		matr_init_rand(10,1000);
	}
	MPI_Bcast(&vert_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
	comp = new int[vert_num];
	sendcounts = new int[ProcNum];
	displs = new int[ProcNum];
	if (!ProcRank)
	{
		mst = new int[vert_num*vert_num];
		for (int i = 0; i < vert_num*vert_num; i++)
			mst[i] = 0;

		for (int i = 0; i < vert_num; i++)
			comp[i] = i;

		

		int vert_proc = vert_num / (ProcNum - 1);
		int vert_add = vert_num % (ProcNum - 1);
		/*if (vert_add)
			sendcounts[0] = (vert_proc + 1) * vert_num;
		else
			sendcounts[0] = vert_proc * vert_num;*/
		sendcounts[0] = 0;
		displs[0] = 0;
		for (int i = 1; i <= vert_add; i++)
		{
			sendcounts[i] = (vert_proc + 1) * vert_num;
			displs[i] = sendcounts[i - 1] + displs[i - 1];
		}
		for (int i = vert_add + 1; i < ProcNum; i++)
		{
			sendcounts[i] = vert_proc * vert_num;
			displs[i] = sendcounts[i - 1] + displs[i - 1];
		}
		/*sendcounts[ProcNum - 1] = 0;
		displs[ProcNum - 1] = vert_num*vert_num;*/
	}
}

void calculate()
{
	for (int i = 0; i < proc_vert_num; i++)
	{
		int min = _CRT_INT_MAX, en_vert = -1;
		for (int j = 0; j < vert_num; j++)
		{
			if (proc_matr[i*vert_num+j] && comp[i+proc_displs] != comp[j] && proc_matr[i*vert_num+j] < min)
			{
				min = proc_matr[i*vert_num+j];
				en_vert = j;
			}
		}
		proc_calc_en_vert[i] = en_vert;
	}
}

bool components_update()
{
	bool update = false;
	os << "Components:" << endl;
	for (int i = 0; i < vert_num; i++)
	{
		os << comp[i] << " ";
	}
	os << endl << "Calculated vertexes:" << endl;
	for (int i = 0; i < vert_num; i++)
	{
		os << calc_en_vert[i] << " ";
	}
	os << endl << "Current MST:" << endl;
	for (int i = 0; i < vert_num; i++)
	{
		for (int j = 0; j < vert_num; j++)
			os << mst[i*vert_num + j] << " ";
		os << endl;
	}
	os << "-----------------------" << endl;
	for (int i = 0; i < vert_num; i++)
	{
		if (calc_en_vert[i] != -1 && comp[i] != comp[calc_en_vert[i]])
		{
			update = true;
			int min = matr[i*vert_num+calc_en_vert[i]], st_vert = i, en_vert = calc_en_vert[i];
			for (int j = i + 1; j < vert_num; j++)
			{
				if (comp[j] == comp[i] && calc_en_vert[j] != -1)
				{
					if (matr[j*vert_num + calc_en_vert[j]] < matr[st_vert*vert_num+en_vert])
					{
						min = matr[j*vert_num + calc_en_vert[j]];
						st_vert = j;
						en_vert = calc_en_vert[j];
					}
					calc_en_vert[j] = -1;
				}
			}
			mst[st_vert*vert_num+en_vert] = mst[en_vert*vert_num+st_vert] = min;

			if (comp[st_vert] < comp[en_vert])
			{
				for (int j = 0; j < vert_num; j++)
					if (comp[j] == comp[en_vert] && j != en_vert)
						comp[j] = comp[st_vert];
				comp[en_vert] = comp[st_vert];
				comp_num--;
				edges_num++;
				f_opt += min;
			}
			else if (comp[st_vert] > comp[en_vert])
			{
				for (int j = 0; j < vert_num; j++)
					if (comp[j] == comp[st_vert] && j != st_vert)
						comp[j] = comp[en_vert];
				comp[st_vert] = comp[en_vert];
				comp_num--;
				edges_num++;
				f_opt += min;
			}
		}
	}
	return update;
}

void opt_count()
{
	f_opt = 0;
	for (int i = 0; i < vert_num; i++)
		for (int j = i + 1; j < vert_num; j++)
			if (matr[i*vert_num + j])
				f_opt += matr[i*vert_num + j];
}

void results_record()
{
	//os.open("results.txt");
	for (int i = 0; i < vert_num; i++)
	{
		for (int j = 0; j < vert_num; j++)
			os << mst[i*vert_num+j] << " ";
		os << endl;
	}
	os << "------------------" << endl;
	os << "Edges number: " << edges_num << endl;
	os << "Components number: " << comp_num << endl;
	if (edges_num == vert_num - comp_num)
		os << "Calculated graph is spanning tree (forest)" << endl;
	else
		os << "Calculated graph is NOT spanning tree (forest). Results are not correct!" << endl;
	//opt_count();
	os << endl << "MST edges sum: " << f_opt << endl;
	os << "Parallel version time: " << en_time - st_time << endl;
	os.close();
}

void mem_del()
{
	delete[] comp;
	delete[] calc_en_vert;

	if (!ProcRank)
	{
		delete[] matr;
		delete[] mst;
		delete[] sendcounts;
		delete[] displs;
	}
}

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if(!ProcRank)
	{
		os.open("results.txt");
	}
	mem_init();
	MPI_Bcast(sendcounts, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	/*MPI_Bcast(displs, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	proc_vert_num = sendcounts[ProcRank];
	proc_displs = displs[ProcRank];*/
	/*MPI_Scatter(sendcounts, ProcNum, MPI_INT, &proc_vert_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(displs, ProcNum, MPI_INT, &proc_displs, 1, MPI_INT, 0, MPI_COMM_WORLD);*/
	proc_matr = new int[sendcounts[ProcRank]];
	MPI_Scatterv(matr, sendcounts, displs, MPI_INT, proc_matr, sendcounts[ProcRank], MPI_INT, 0, MPI_COMM_WORLD);

	if (!ProcRank)
	{
		for (int i = 0; i < ProcNum; i++)
		{
			sendcounts[i] /= vert_num;
			displs[i] /= vert_num;
		}
		comp_num = vert_num;
		edges_num = 0;
		f_opt = 0;
	}

	MPI_Bcast(sendcounts, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displs, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	proc_vert_num = sendcounts[ProcRank];
	proc_displs = displs[ProcRank];

	if (ProcRank)
		proc_calc_en_vert = new int[proc_vert_num];
	else
		calc_en_vert = new int[vert_num];

	MPI_Barrier(MPI_COMM_WORLD);
	st_time = MPI_Wtime();
	while (true)
	{
		MPI_Bcast(comp, vert_num, MPI_INT, 0, MPI_COMM_WORLD);

		if (ProcRank)
		{
			calculate();
		}
		MPI_Gatherv(proc_calc_en_vert, proc_vert_num, MPI_INT, calc_en_vert, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
		
		int update;
		if (!ProcRank)
		{
			update = components_update();
		}
		MPI_Bcast(&update, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (!update)
		{
			os << "Cycle exit" << endl;
			break;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	en_time = MPI_Wtime();
	if (!ProcRank)
	{
		os << "Cycle exit" << endl;
	}
	if (!ProcRank)
	{
		os << "Results record start" << endl;
		results_record();
	}

	mem_del();
	MPI_Finalize();
	return 0;
}