#include <stdio.h>
#include <fstream>
#include <stack>
using namespace std;

int **matr;
int **mst;
int* comp_num;
int vert_num, components, f_opt;
stack<int> curr_edges;

void init_file(char* file_name)
{
	ifstream is;
	is.open(file_name);
	is >> vert_num;
	matr = new int*[vert_num];
	for (int i = 0; i < vert_num; i++)
		matr[i] = new int[vert_num];

	for (int i = 0; i < vert_num; i++)
	{
		for (int j = 0; j < vert_num; j++)
		{
			is >> matr[i][j];
			//printf("%d ", matr[i][j]);
		}
		//printf("\n");
	}
	is.close();
}

void calculate()
{
	f_opt = 0;
	//components = vert_num;
	while (true)
	{
		for (int i = 0; i < vert_num; i++)
		{
			int min = _CRT_INT_MAX, st_vert = -1, en_vert = -1;
			for (int j = 0; j < vert_num; j++)
			{
				if (comp_num[j] == i)
				{
					for (int k = 0; k < vert_num; k++)
					{
						if (matr[k][j] && comp_num[k] != comp_num[j] && matr[k][j] < min)
						{
							min = matr[k][j];
							st_vert = j;
							en_vert = k;
						}
					}
				}
			}
			if (en_vert != -1)
			{
				curr_edges.push(en_vert);
				curr_edges.push(st_vert);
			}
		}
		if (curr_edges.empty())
			break;
		while (!curr_edges.empty())
		{
			int st_vert = curr_edges.top();
			curr_edges.pop();
			int en_vert = curr_edges.top();
			curr_edges.pop();

			if (mst[st_vert][en_vert] != matr[st_vert][en_vert])
			{
				mst[st_vert][en_vert] = mst[en_vert][st_vert] = matr[st_vert][en_vert];
				//components--;
				f_opt += mst[st_vert][en_vert];
				int curr_comp, new_comp;
				if (comp_num[st_vert] > comp_num[en_vert])
				{
					curr_comp = comp_num[st_vert];
					new_comp = comp_num[en_vert];
				}
				else
				{
					curr_comp = comp_num[en_vert];
					new_comp = comp_num[st_vert];
				}
				for (int i = 0; i < vert_num; i++)
					if (comp_num[i] == curr_comp)
						comp_num[i] = new_comp;
			}
		}
	}
}

void results_record(char* file_name)
{
	ofstream os;
	os.open(file_name);
	for (int i = 0; i < vert_num; i++)
	{
		for (int j = 0; j < vert_num; j++)
			os << mst[i][j] << " ";
		os << endl;
	}
	os << "------------------" << endl;
	os << "MST edges sum: " << f_opt << endl;
	os.close();
}

void mem_init()
{
	mst = new int*[vert_num];
	for (int i = 0; i < vert_num; i++)
		mst[i] = new int[vert_num];

	for (int i = 0; i < vert_num; i++)
		for (int j = 0; j < vert_num; j++)
			mst[i][j] = 0;

	comp_num = new int[vert_num];
	for (int i = 0; i < vert_num; i++)
		comp_num[i] = i;
}

void mem_del()
{
	for (int i = 0; i < vert_num; i++)
		delete[] mst[i];
	delete[] mst;
	for (int i = 0; i < vert_num; i++)
		delete[] matr[i];
	delete[] matr;
	delete[] comp_num;
}

int main(int argc, char** argv)
{
	init_file("adj_matr.txt");
	mem_init();

	calculate();

	results_record("results.txt");
	mem_del();

	return 0;
}