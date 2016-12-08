#include <stdio.h>
#include "mpi.h"
#include <fstream>
#include <time.h>
#include <stdlib.h>
using namespace std;

int *matr=nullptr;
int *mst_seq=nullptr, *mst_par=nullptr;
int *comp=nullptr;
int *calc_en_vert=nullptr;
int vert_num, comp_num;
int f_opt_seq, f_opt_par;
int edges_num;
ofstream os;
ifstream is;

int ProcNum, ProcRank;
int *sendcounts=nullptr, *displs=nullptr;
int *proc_matr=nullptr;
int *proc_calc_en_vert=nullptr;
int proc_vert_num, proc_displs;

double st_time_seq, en_time_seq;
double st_time_par, en_time_par;
MPI_Status status;

void matr_init()                             //заполнение матрицы из файла
{
	is.open("adj_matr.txt");
	is >> vert_num;
	matr = new int[vert_num*vert_num];

	for (int i = 0; i < vert_num; i++)
		for (int j = 0; j < vert_num; j++)
			is >> matr[i*vert_num + j];

	is.close();
}

void matr_init_rand(int v_num, int power)      //заполнение матрицы случайными числами
{
	srand(time(nullptr));
	vert_num = v_num;
	matr = new int[vert_num*vert_num];
	for (int i = 0; i < vert_num; i++)
	{
		matr[i*vert_num + i] = 0;
		for (int j = i + 1; j < vert_num; j++)
		{
			matr[i*vert_num + j] = matr[j*vert_num + i] = rand() % power;
		}
	}

	if (vert_num < 30)
	{
		os << "Adjacency matrix:" << endl;
		for (int i = 0; i < vert_num; i++)
		{
			for (int j = 0; j < vert_num; j++)
				os << matr[i*vert_num + j] << " ";
			os << endl;
		}
		os << "--------------------" << endl;
	}
}

void mem_init()                                         //инициализация памяти
{
	if (!ProcRank)
	{
		matr_init_rand(8000, 1000);
	}
	MPI_Bcast(&vert_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
	comp = new int[vert_num];
	sendcounts = new int[ProcNum];
	displs = new int[ProcNum];
	if (!ProcRank)
	{
		mst_seq = new int[vert_num*vert_num];
		mst_par = new int[vert_num*vert_num];
		calc_en_vert = new int[vert_num];

		for (int i = 0; i < vert_num*vert_num; i++)
			mst_seq[i] = mst_par[i] = 0;

		for (int i = 0; i < vert_num; i++)
			comp[i] = i;

		int vert_proc = vert_num / (ProcNum - 1);
		int vert_add = vert_num % (ProcNum - 1);
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
	}
}



void calculate()                                         //вычисление минимального ребра для каждой вершины в процессе
{
	for (int i = 0; i < proc_vert_num; i++)
	{
		int min = _CRT_INT_MAX, en_vert = -1;
		for (int j = 0; j < vert_num; j++)
		{
			if (proc_matr[i*vert_num + j] && comp[i + proc_displs] != comp[j] && proc_matr[i*vert_num + j] < min)
			{
				min = proc_matr[i*vert_num + j];
				en_vert = j;
			}
		}
		proc_calc_en_vert[i] = en_vert;
	}
}

bool components_update(int* mst,int &f_opt)             //вычисление минимального ребра для каждой компоненты - пересчёт принадлежности вершин компонентам
{
	bool update = false;
	for (int i = 0; i < vert_num; i++)
	{
		if (calc_en_vert[i] != -1 && comp[i] != comp[calc_en_vert[i]])
		{
			update = true;
			int min = matr[i*vert_num + calc_en_vert[i]], st_vert = i, en_vert = calc_en_vert[i];
			for (int j = i + 1; j < vert_num; j++)
			{
				if (comp[j] == comp[i] && calc_en_vert[j] != -1)
				{
					if (matr[j*vert_num + calc_en_vert[j]] < matr[st_vert*vert_num + en_vert])
					{
						min = matr[j*vert_num + calc_en_vert[j]];
						st_vert = j;
						en_vert = calc_en_vert[j];
					}
					calc_en_vert[j] = -1;
				}
			}
			mst[st_vert*vert_num + en_vert] = mst[en_vert*vert_num + st_vert] = min;

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

void sequential_version()                          //последовательная версия
{
	while (true)
	{
		for (int i = 0; i < vert_num; i++)
		{
			int min = _CRT_INT_MAX, en_vert = -1;
			for (int j = 0; j < vert_num; j++)
			{
				if (matr[i*vert_num + j] && comp[i] != comp[j] && matr[i*vert_num + j] < min)
				{
					min = matr[i*vert_num + j];
					en_vert = j;
				}
			}
			calc_en_vert[i] = en_vert;
		}
		if (!components_update(mst_seq, f_opt_seq))
			break;
	}
}

void opt_count(int* mst, int &f_opt)
{
	f_opt = 0;
	for (int i = 0; i < vert_num; i++)
		for (int j = i + 1; j < vert_num; j++)
			if (mst[i*vert_num + j])
				f_opt += mst[i*vert_num + j];
}

bool check()                           //проверка эквивалентности полученных результатов
{
	/*if (f_opt_par != f_opt_seq)
		return false;*/
	for (int i = 0; i < vert_num; i++)
		for (int j = i + 1; j < vert_num; j++)
			if (mst_seq[i*vert_num + j] != mst_par[i*vert_num + j])
			{
				os << "Results are not equal in (i,j)" << endl;
				return false;
			}
	return true;
}

void results_record()                 //запись результатов
{
	if (vert_num < 30)
	{
		os << "Sequential version MST:" << endl;
		for (int i = 0; i < vert_num; i++)
		{
			for (int j = 0; j < vert_num; j++)
				os << mst_seq[i*vert_num + j] << " ";
			os << endl;
		}
		os << "-----------------------" << endl;
		os << "Parallel version MST:" << endl;
		for (int i = 0; i < vert_num; i++)
		{
			for (int j = 0; j < vert_num; j++)
				os << mst_par[i*vert_num + j] << " ";
			os << endl;
		}
	}

	//opt_count(mst_seq, f_opt_seq);
	os << "MST edges sum sequential: " << f_opt_seq << endl;
	os << "MST edges sum parallel: " << f_opt_par << endl;
	if (check())
		os << "Results are equal" << endl;
	else
		os << "Results are not equal" << endl;
	os << "Edges number: " << edges_num << endl;
	os << "Components number: " << comp_num << endl;
	if (edges_num == vert_num - comp_num)
		os << "Calculated graph is spanning tree (forest)" << endl;
	else
		os << "Calculated graph is NOT spanning tree (forest). Results are not correct!" << endl;
	//opt_count(mst_par, f_opt_par);
	os << endl << "Sequential version time: " << en_time_seq - st_time_seq << endl;
	os << "Parallel version time: " << en_time_par - st_time_par << endl;

	double spdup = (en_time_seq - st_time_seq) / (en_time_par - st_time_par);
	os << "Speedup: " << spdup << endl;
	os << "Efficiency: " << spdup / ProcNum << endl;
	os.close();
}

void mem_del()               //освобождение памяти
{
	delete[] proc_matr;
	delete[] comp;
	delete[] calc_en_vert;
	delete[] proc_calc_en_vert;

	if (!ProcRank)
	{
		delete[] matr;
		delete[] mst_seq;
		delete[] mst_par;
		delete[] sendcounts;
		delete[] displs;
	}
}

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (!ProcRank)
	{
		os.open("results.txt");
	}

	mem_init();

	if (!ProcRank)                       //последовательная версия
	{
		f_opt_seq = 0;
		st_time_seq = MPI_Wtime();
		sequential_version();
		en_time_seq = MPI_Wtime();               //конец последовательной версии

		for (int i = 0; i < vert_num; i++)
			comp[i] = i;
	}

	MPI_Bcast(sendcounts, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	proc_matr = new int[sendcounts[ProcRank]];                                    //распределение вершин между процессами
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
		f_opt_par = 0;
	}

	MPI_Bcast(sendcounts, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displs, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	proc_vert_num = sendcounts[ProcRank];
	proc_displs = displs[ProcRank];

	if (ProcRank)
		proc_calc_en_vert = new int[proc_vert_num];

	MPI_Barrier(MPI_COMM_WORLD);
	st_time_par = MPI_Wtime();
	while (true)                       //параллельная версия
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
			update = components_update(mst_par, f_opt_par);
		}
		MPI_Bcast(&update, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (!update)
			break;
	}                                     //конец параллельной версии
	MPI_Barrier(MPI_COMM_WORLD); 
	en_time_par = MPI_Wtime();

	if (!ProcRank)
	{
		results_record();
	}

	mem_del();
	MPI_Finalize();
	return 0;
}