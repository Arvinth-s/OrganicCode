//OM NAMO NARAYANA
#include<iostream>
#include<stdio.h>
#include<list>
#include<vector>
#include<string>
#include<math.h>
#include<limits>

/*
	Concepts of
		1)trees
		2)graphs
		3)searching
		4)sorting
		5)hashing
	has been implemented in this algorithm

	The input is a unidirectional map which points from reaction to product at each step
	The program predicts the yeild of each and every products/intermediates
	The input should also have the information about the order of the reaction
	We directly provide the integrated form of the reaction rate and calculate the concentrations
	Constraints:The programs expects the users to feed in the data for energy of the intermediates 
	possible solution: We can build a file or hash table that stores the new data and avoids the user to input the energy of common intermediates
*/

/*
	The user has multiple choices
		1)He can directly provide the value of Kf and Kb
		2)He can input the Temperature and Gibbs free energy of the compounds involved in the reaction
		3)He can input the Tempearture and change in  Gibbs free energy of the compounds
		4)He can input the totalEnergy(H) change, entropy change and temperature and thus change in Gibbs free energy can be calculated
		
*/

using namespace std;

#define R 8.31446261815324
#define unassigned 0.0012317

struct hashTable{
	string key;
	double energy;
	double state;
	char state;
	string allotropy;
};

int findKeyforHash(string str)
{
	
}

struct Compound{
	string name;
	double energy;
	double entropy;
	list<struct Compoud> *products;
	list<struct Compound> *reactants;
	//list<struct compound> *friends;
};

compound::compound(string name, double energy = unassigned, double entropy = unassigned)
{
	this->name = name;
	this->energy = energy;
	this->entropy = entropy;
	products = NULL;
	reactants = NULL;
	//friends = NULL;
}

struct Reaction{
	Compound **
};

double KwithG(double dG, float T)
	return (exp(-dG /(R*T) ));
	
double GwithHS(double dH, double dS, float T)
	return (dH - T*dS);

double CalConc(float n, double conc, float time, float k)
{
	if(n == 1)
		return (exp(-1*k*time)*conc);
	else if(n == 0)
		return (conc - k*time);
	else
	{
		double step1 = k*time + 1/pow(conc, n-1);
		return (1/pow(step1, 1/(n-1)));
	}
}

double findK(float n, float t0, double conc0, float t1, double conc1)
{
	if(n == 1)
		return(log(conc1/conc2)/(t2 - t1));
	else if(n==0)
		return (conc1 - conc2)/(t2-t1);
	else
		return (1/pow(conc2, n-1) - 1/pow(conc1, n-1))/t;
	
}

double loss(list<double> predicted, list<double> conc)
{
	double sum = 0;
	list<double>::iterator itr;
	list<double>::iterator jtr = conc.begin();	
	for(itr = predicted.begin(); itr!=predicted.end(); itr++, jtr++)
		sum += pow(*itr - *jtr, 2);
	return sum;
}

//This function returns the order of the reaction with precision 0.5 given the concentration of the reactants at more than 2 time
float calcOrder(list<float> time, list<double> conc)
{
	if(time.size() != conc.size())
		printf("error:mismatch in size\n"),return -1.0;
	if(time.size() <= 2)
		printf("error:inadequate data, provide atleast 3 data\n"), return -1.0;
	list<float>::iterator titr;
	list<double>::iterator citrl
	double loss = numeric_limits<T>::max();
	double closs = 0;
	float count = 0;
	titr = time.begin(), citr = conc.begin();
	t0 = *titr, conc0 = *citr;
	titr++, citr++;
	t1 = *titr, conc1 = *citr;
	int size = time.size();
	citr--;
	while(loss >= closs)
	{
		loss = closs;
		double k = findK(count, t0, conc0, t1, conc1);
		list<double> predicted;
		predicted.push_back(*citr);
		for(int i = 0; i < size; i++)
			predicted.push_back(CalConc(count, predicted.last(), *titr, k)), titr++; 
		closs = loss(predicted, conc);
		count += 0.5;
	}
	return count;
}

int main()
{
	float order;
	printf("Enter the order of the reaction:");
	scanf("%f", &order);
	double conc;
	float time, k;
	printf("Enter the initial concentration of the reaction:");
	scanf("%lf", &conc);
	printf("Enter the inital rate constant of the reaction:");
	scanf("%f", &k);
	printf("Enter the time at which you wish to find the concentration of the reactant:");
	scanf("%f", &time);
	double pro = CalConc(order, conc, time, k);
	cout<<pro;
	return 0;
}
