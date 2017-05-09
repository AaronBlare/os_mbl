#define _CRT_SECURE_NO_DEPRECATE
#include "header.h"


void read_matrix_bin (FILE * file, MKL_Complex16 * matrix, unsigned int N)
{
	fread(matrix, sizeof(MKL_Complex16), N * N, file);
}

void read_branch_bin (FILE * file, split * branch)
{
	unsigned int deep, N;
	split * node = branch;
	fread(&deep, sizeof(int), 1, file);
	N = (branch->prev)->N;
	for(unsigned int i = 0; i < deep - 1; i++)
	{
		fread(&(node->dt), sizeof(double), 1, file);
		fread(&(node->steps), sizeof(int), 1, file);
		node->counter = 0;
		node->g = 0;
		node->N = N;
		node->type = true;
		node->matrix = new MKL_Complex16[N * N];
		read_matrix_bin(file, node->matrix, N);
		node->next = new split;
		node->next->prev = node;
		node = node->next;
	}
	fread(&(node->dt), sizeof(double), 1, file);
	fread(&(node->steps), sizeof(int), 1, file);
	node->counter = 0;
	node->g = 0;
	node->N = N;
	node->type = true;
	node->matrix = new MKL_Complex16[N * N];
	read_matrix_bin(file, node->matrix, N);
	node->next = 0;
}

split * create_struct_bin (FILE * file)
{
	split * head = new split[1];
	unsigned int q_branch, N;
	double T;
	fread(&(T), sizeof(double), 1, file);
	fread(&(N), sizeof(int), 1, file);
	fread(&(q_branch), sizeof(int), 1, file);
	head->prev = 0;
	head->type = false;
	head->dt = T;
	head->counter = q_branch;
	head->N = N;
	head->next =  new split[q_branch];
	for(unsigned int i = 0; i < q_branch; i++)
	{
		(head->next)[i].prev = head;
		read_branch_bin(file, &((head->next)[i]));
	}

	fread(&(head->steps), sizeof(int), 1, file);;
	head->matrix = new MKL_Complex16[head->steps * N * N];
	head->g = new double[head->steps];
	for(unsigned int k = 0; k < head->steps; k++)
	{
		fread(&(head->g[k]), sizeof(double), 1, file);
		read_matrix_bin(file, &(head->matrix[k * N * N]), N);
	}
	return head;
}


void read_matrix (FILE * file, MKL_Complex16 * matrix, unsigned int N)
{
	for(unsigned int i = 0; i < N; i++)
	{
		for(unsigned int j = 0; j < N; j++)
		{
			fscanf(file, "%lf %lf", &matrix[i * N + j].real, &matrix[i * N + j].imag);
		}
	}
}

void read_branch (FILE * file, split * branch)
{
	unsigned int deep, N;
	split * node = branch;
	fscanf (file, "%d", &deep);
	N = (branch->prev)->N;
	for(unsigned int i = 0; i < deep - 1; i++)
	{
		fscanf(file, "%lf %d", &(node->dt), &(node->steps));
		node->counter = 0;
		node->g = 0;
		node->N = N;
		node->type = true;
		node->matrix = new MKL_Complex16[N * N];
		read_matrix(file, node->matrix, N);
		node->next = new split;
		node->next->prev = node;
		node = node->next;
	}
	fscanf(file, "%lf %d", &(node->dt), &(node->steps));
	node->counter = 0;
	node->g = 0;
	node->N = N;
	node->type = true;
	node->matrix = new MKL_Complex16[N * N];
	read_matrix(file, node->matrix, N);
	node->next = 0;
}

split * create_struct (FILE * file)
{
	split * head = new split[1];
	unsigned int q_branch, N;
	double T;
	fscanf (file, "%lf", &T);
	fscanf (file, "%d", &N);
	fscanf (file, "%d", &q_branch);
	head->prev = 0;
	head->type = false;
	head->dt = T;
	head->counter = q_branch;
	head->N = N;
	head->next =  new split[q_branch];
	for(unsigned int i = 0; i < q_branch; i++)
	{
		(head->next)[i].prev = head;
		read_branch(file, &((head->next)[i]));
	}

	fscanf (file, "%d", &head->steps);
	head->matrix = new MKL_Complex16[head->steps * N * N];
	head->g = new double[head->steps];
	for(unsigned int k = 0; k < head->steps; k++)
	{
		fscanf (file, "%lf", &(head->g[k]));
		read_matrix(file, &(head->matrix[k * N * N]), N);
	}
	return head;
}

void cmp_branch_not_member(split * src, split * dst)
{
	split * node1 = src;
	split * node2 = dst;
	while (node1->next)
	{
		node2->dt = node1->dt;
		node2->steps = node1->steps;
		node2->counter = node1->counter;
		node2->g = node1->g;
		node2->N = node1->N;
		node2->type = node1->type;
		node2->matrix = node1->matrix;
		node2->next = new split;
		node2->next->prev = node2;
		node2 = node2->next;
		node1 = node1->next;
	}
	node2->dt = node1->dt;
	node2->steps = node1->steps;
	node2->counter = node1->counter;
	node2->g = node1->g;
	node2->N = node1->N;
	node2->type = node1->type;
	node2->matrix = node1->matrix;
	node2->next = 0;
}

void cmp_struct_not_member(split * src, split * dst)
{
	dst->prev = src->prev;
	dst->type = src->type;
	dst->dt = src->dt;
	dst->counter = src->counter;
	dst->N = src->N;
	dst->next =  new split[dst->counter];
	for(unsigned int i = 0; i < dst->counter; i++)
	{
		(dst->next)[i].prev = dst;
		cmp_branch_not_member(&((src->next)[i]), &((dst->next)[i]));
	}

	dst->steps = src->steps;
	dst->matrix = src->matrix;
	dst->g = src->g;
}

void delete_branch (split * branch)
{
	delete (branch->matrix);
	branch->prev = 0;
	branch->matrix = 0;
	branch->steps = 0;
	branch->type = false;
	branch->dt = 0;
	branch->counter = 0;
	branch->N = 0;

	branch = branch->next;
	while (branch->next->next == 0)
	{
		delete (branch->matrix);
		branch->prev = 0;
		branch->matrix = 0;
		branch->steps = 0;
		branch->type = false;
		branch->dt = 0;
		branch->counter = 0;
		branch->N = 0;

		branch = branch->next;
		delete (branch->prev);
	}
	delete (branch->matrix);
	branch->prev = 0;
	branch->matrix = 0;
	branch->steps = 0;
	branch->type = false;
	branch->dt = 0;
	branch->counter = 0;
	branch->N = 0;

	delete (branch);
}

void delete_split_struct (split * head)
{
	for(unsigned int i = 0; i < head->counter; i++)
	{
		delete_branch (&(head->next)[i]);
	}
	delete (head->next);
	delete (head->matrix);
	delete (head->g);
}

void delete_branch_not_member (split * branch)
{
	branch->prev = 0;
	branch->matrix = 0;
	branch->steps = 0;
	branch->type = false;
	branch->dt = 0;
	branch->counter = 0;
	branch->N = 0;

	branch = branch->next;
	while (branch->next->next == 0)
	{
		branch->prev = 0;
		branch->matrix = 0;
		branch->steps = 0;
		branch->type = false;
		branch->dt = 0;
		branch->counter = 0;
		branch->N = 0;

		branch = branch->next;
		delete (branch->prev);
	}
	branch->prev = 0;
	branch->matrix = 0;
	branch->steps = 0;
	branch->type = false;
	branch->dt = 0;
	branch->counter = 0;
	branch->N = 0;

	delete (branch);
}

void delete_split_struct_not_member (split * head)
{
	for(unsigned int i = 0; i < head->counter; i++)
	{
		delete_branch_not_member (&(head->next)[i]);
	}
	delete (head->next);
	head->matrix = 0;
	head->g = 0;
}