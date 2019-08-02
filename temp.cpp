//OM NAMO NARAYANA
#include<iostream>
#include<stdio.h>
#include<vector>
#include<list>
#include<string.h>
#include<string>


//It All elements having more than one valency can form more than one bond
//So we need to find a notation  to  represent every atom uniquely to add elements to that atom
//For this case we use a hashtable


using namespace std;
char h1[] = "H1";
int count = -1;


int encode(char ele[], int size);
void QuickSort(int *arr, int *arange, int start, int end);
vector<int*> SetFlags(int *arr, int *arange);
int Search(int *arr, int ele, vector<int*> flags);
int GetHash(int code);
void PrintFlags(vector<int*> flags);
void PrintArray(int *arr, int n);
void PrintVector(vector<int> arr);
string Decode(int code);

struct compound{
	static int count;
	int key;//index of the object of the structure corresponding to this atom
	char atom[4];//Chemical formula for that atom
	int V;//Maximum possible valency of that atom
	int state;//Oxidation state
	int code;//integer identity
	vector<int> element;//Contains the list of bonded atoms' code
	vector<int> mbonds;
	int electronegativity;//Electronegativity of that atom
	int filled;//This is used to count the number of atoms filled
	int charge;
	int lonepair;
	int Bond(char batom[]);//To form a bond
	void BreakBond(int code, int flag);//To break a bond
	void CreateAtom(char atom[]);//constructor
	int CalcElectroNegativity(int code);
	//int CalcElectroNegativity();
};

//This will return the index of the least number that is greater than or equal to that number and whether that series is continuous and the number of elements in that sequence
vector<int> Seq(int *arr, int ele, int start)
{
	int i = start;
	int count = 1;
	int flag = 0;
	while(int(arr[start]/100) == int(arr[i]/ 100))
	{
		if(arr[i]%100 == count)
			count++;
		else {
		flag = -1;	}
		i++;
	}
	vector<int> temp;
	temp.push_back(flag);
	temp.push_back(i - start - 1);
	temp.push_back(ele);
	return temp;
}

void PrintVector(vector<int> arr)
{
	vector<int>::iterator itr;
	cout<<"Printing vector\n";
	for(itr = arr.begin(); itr!= arr.end(); itr++)
		cout<<*itr<<"\t";
	cout<<endl;
}

vector<vector<int> > SetFlags(int *arr, int *arange, int n)
{
	int start = 0;
	QuickSort(arr, arange, start, n-1);
	vector<int*> flag;
	vector<int> temp_vector;
	vector<vector<int> > index;
	for(int i = 0; i < n; i++)
	{
		index.push_back(Seq(arr, arr[i], i));
		temp_vector = index.back();
		i += temp_vector[1];
	}
	return index;
}


void PrintArray(int *arr, int n)
{
	printf("Printing the array\n");
	for(int i = 0; i < n; i++)
		cout<<arr[i]<<"\t";
	printf("\n");
	return;
}

void PrintFlags(vector<vector<int> > flags)
{
	cout<<"Printing flags\n";
	vector<vector<int> >::iterator itr;
	vector<int> temp;
	for(itr = flags.begin(); itr!= flags.end(); itr++)
	{
		temp = *itr;
		printf("%d\t%d\t%d\n", temp[0], temp[1], temp[2]);
	}
}

void swap(int *arr, int i, int j)
{
	int temp = arr[i];
	arr[i] = arr[j];
	arr[j] = temp;
}

//end is an included index
void QuickSort(int *arr, int *arange, int start, int end)
{
	if(start >= end)
		return;
	int i = start - 1, j = start;
	for(j = start; j < end; j++)
		if(arr[j] < arr[end])
		{
			i++;
			swap(arr, i, j);
			swap(arange, i , j);
		}
	i++;
	if(i < end)
	{
		swap(arr, i, end);
		swap(arange, i, end);
	}
	if(i-1 > start)QuickSort(arr, arange, start, i-1);
	if(end > i+1)QuickSort(arr, arange, i+1, end);
	return;
}

int Search(int *arr, int ele, vector<vector<int> > flags)
{
	int sum = 0;
	vector<vector<int> >::iterator itr;
	vector<int> temp_vector;
	for(itr = flags.begin(); itr!=flags.end();itr++)
	{
		temp_vector = *itr;
		if(temp_vector[2]/100 == ele/100)
		{
			if(!temp_vector[0])
			{
				return (sum + ele%100);
			}
			else
			{
				for(int k = 0; k< temp_vector[1]+1; k++)
				{
					if(arr[sum+k] == ele)
					{
						return sum+k;
					}
				}
			}
		}
		else if(temp_vector[2]/100 > ele/100)
		{
			itr--;
			temp_vector = *itr;
			sum -= temp_vector[1]+1;
			if(!temp_vector[0])
			{
				return (sum + ele%100);
			}
			else
			{
				for(int k = 0; k< temp_vector[1]+1; k++)
				{
					if(arr[sum+k] == ele)
					{
						return sum+k;
					}
				}
			}
		}
		sum += temp_vector[1]+1;//index of the next element
	}
	cout<<"leaving search with negative 1\n";
	return -1;
}

int GetHash(int *arr, int *arange, int code, vector<vector<int> > flags, int n)
	{
		int temp_hold = Search(arr, code, flags);
		if(!temp_hold)
		{
			printf("Something went wrong. Please check the input. Terminating the Program\n");
			return -1;
		}
		else
		{
			return temp_hold;
		}
	}

//This function is supposed to encode the string to int for a compound identity
int encode(char ele[])
{
	int i=0;
	char element[3];
	int pos = 0;
	if(strlen(ele) == 2)
	{
		pos = ele[1] - '0'; 
		element[0] = ele[0];
		element[1] = '\0';
	}
	else if(strlen(ele) == 2)
	{
		pos = ele[2] - '0'; 
		element[0] = ele[0];
		element[1] =  ele[1] ;
	}
	else
	{
		printf("Invalid choice, The atom should be of form XXY, where XX should be character and Y should be a natural number\n"); 
		return -1;
	}
	
	char g1[][3] = {"H", "Li", "Na", "K", "Rb"};
	char halogens[][3] = {"F", "Cl", "Br", "I"};
	char NFam[][3] = {"N", "P", "As", "Sb", "Bi"};
	char OFam[][3] = {"O", "S", "Se", "Te", "Po"};
	char BFam[][3] = {"B", "Al"};
	char carbon[][3] = {"C"};
	for(i = 0; i < 5; i++)
		if(!strcmpi(element, g1[i]))
			return 1000+i*100+pos;
	for(i = 0; i < 4; i++)
		if(!strcmpi(element, halogens[i]))
			return 7000+100*i+pos;
	for(i = 0; i < 5; i++)
		if(!strcmpi(element, NFam[i]))
			return 5000+i*100+pos;
	for(i = 0; i < 5; i++)
		if(!strcmpi(element, OFam[i]))
			return 6000+i*100+pos;
	for(i = 0; i < 2; i++)
		if(!strcmpi(element, BFam[i]))
			return 3000+i*100+pos;
	if(!strcmpi(element, carbon[0]))
		return 4000+pos;
	
}

string Decode(int code)
{
	char group[][5][3] = {{"B ", "Al", "Ga", "In", "Tl"}, {"C ", "Si", "Ge", "Sn","Pb"}, {"N ", "P ", "As", "Sb", "Bi"}, {"O ", "S ", "Se", "Te", "Po"}, {"F ", "Cl", "Br", "I ", "At"}};
	char Hgroup[][3] = {"H", "Li", "Na", "K", "Rb"};
	char tempH[] = " H";
	string atom = "";
	if(code/1000 >1)
	{
		atom = atom + group[code/1000 - 3][(code/100)%10][0];
		atom = atom + group[code/1000 - 3][(code/ 100)%10][1];
	}
	else
	{
		
		atom = atom+tempH;
	}
	return atom;
}
//This class contains information of each atom
//Each atom in the compound has its object in compound class, i.e., no of atoms = no. of objects
//V represents the valency of the compounds, since each object represents an atom we can set the maximum no of bonds for that atom



int compound::count = 0;

int compound::CalcElectroNegativity(int co)
{
	if(!co)
		return (this->code)/1000 - 2 *((this->code)%1000-(this->code)%100)/100;
	else
		return (co)/1000 - 2 *((co)%1000-(co)%100)/100;
}

void compound::CreateAtom(char atom[])
{
	strcpy(this->atom, atom);
	this->code = encode(this->atom);
	this->V = (this->code)/1000;
	(V > 5 && (this->code / 100)%10 == 0)?V = 8 - V :V = V;
	this->state = 0;
	this->electronegativity = CalcElectroNegativity(0);
	this->charge = 0;
	this->lonepair = V;
}


void PrintCompound(struct compound *com, int n)
{
	printf("Printing the compound\n");
	for(int i = 0; i<n; i++)
	{
		cout<<"\n\n\t\t\t\t atom:"<<com[i].atom<<endl;
		cout<<"state: "<<com[i].state<<"\t  relative electronegativity:"<<com[i].electronegativity<<"\t  valency:"<<com[i].V<<"\t  charge: "<<com[i].charge<<endl;
		cout<<"It is connected to    atoms 		bonds\n";
		for(int j = 0; j < com[i].element.size(); j++)
			cout<<"                        "<<Decode(com[i].element[j])<<" 		"<<com[i].mbonds[j]<<endl;
		
	}
}


//This adds the bonds between the atoms and updates the electronegativity of the atom
int compound::Bond(char batom[])
{
	if(this->element.size() > V)
	{
		printf("The compound has reached its maximum valency\n");
		return 1;
	}
	for(int i = 0; i < this->element.size(); i++)
	{
		if(encode(batom) == element[i])
		{
			this->mbonds[i] = this->mbonds[i]+1;
			if (this->electronegativity > CalcElectroNegativity(encode(batom)))
				this->state = this->state-1;
			else if(this->electronegativity < CalcElectroNegativity(encode(batom)))
				this->state = this->state + 1;
			else
				printf("Both have same electronegativity\t%d\t%d\t%d\n", this->electronegativity, CalcElectroNegativity(encode(batom)), encode(batom));
			this->lonepair = this->lonepair - 1;
			return 1;
		}
	}
	this->mbonds.push_back(1);
	this->element.push_back(encode(batom));
	if (this->electronegativity > CalcElectroNegativity(encode(batom)))
		this->state = this->state-1;
	else if(this->electronegativity < CalcElectroNegativity(encode(batom)))
		this->state = this->state + 1;
	else
		printf("THIS ATOM: %s \tBoth have same electronegativity\t%d\t%d\t%d\n", this->atom, this->electronegativity, CalcElectroNegativity(encode(batom)), encode(batom));
	return 0;
}

void compound::BreakBond(int code, int flag)
{
	list<int> stack;
	list<int> bonds_stack;
		while(this->element.size() > 0)
		{
			if(this->element.back() == code)
			{
				if(this->mbonds.back() == 1)
				{
					this->element.pop_back();
					this->mbonds.pop_back();
				}
				else
					this->mbonds.back() -= 1;
				if(flag)
				{
					if(this->electronegativity > CalcElectroNegativity(code))
						this->charge -= 1;
					else if(this -> electronegativity == CalcElectroNegativity(code))
						this->charge -= 1;
					else
						this->charge += 1;
				}
				else
				{
					if(this->electronegativity > CalcElectroNegativity(code))
						this->charge += 1;
					else if(this -> electronegativity == CalcElectroNegativity(code))
						this->charge += 1;
					else
						this->charge -= 1;
				}
				break;
			}
			else
			{
				stack.push_back(this->element.back());
				bonds_stack.push_back(this->mbonds.back());
				this->element.pop_back();
				this->mbonds.pop_back();
			}
		}
		while(stack.size() > 0)
		{
			this->element.push_back(stack.back());
			this->mbonds.push_back(bonds_stack.back());
			stack.pop_back();
			bonds_stack.pop_back();
		}
		return;
}

//Iterative BFS is used to calculate resonance
void resonance(struct compound *model, int* nodes, int*arange, vector<vector<int> > flags, int count_bonds, int n_atom, int mbonds)
{
	//cout<<"The number of sigma bonds in the model = "<<count_bonds<<endl;
	//model[0].BreakBond(model[0].element[0], 1);
	//model[arange[GetHash(nodes, arange, model[0].element[0], flags, n_atom) - 1]].BreakBond(model[0].code, 1);
	//(model, n_atom);
	bool visited[n_atom];
	for(int i = 0; i < n_atom; i++)
		visited[i] = false;
	list<struct compound> queue;
	compound current;  
	queue.push_back(model[0]);
	visited[0] = true;
	while(!queue.empty())
	{
		current = queue.front();
		queue.pop_front();
		visited[arange[GetHash(nodes, arange, current.code, flags, n_atom)]] = true;
		int flag = 0;
		for(int i = 0; i < current.mbonds.size(); i++)
		{
			if(current.mbonds[i] >= 2)
			{
				for(int j = 0; j < model[arange[GetHash(nodes, arange, current.element[i], flags, n_atom) - 1]].mbonds.size(); j++)
				for(int k = 0; k < model[arange[GetHash(nodes, arange, model[arange[GetHash(nodes, arange, current.element[i], flags, n_atom) - 1]].element[j], flags, n_atom) - 1]].mbonds.size(); k++)
					if(model[arange[GetHash(nodes, arange, current.element[i], flags, n_atom) - 1]].mbonds[j] >= 2)
					{
						flag = 1;
						model[arange[GetHash(nodes, arange, current.code, flags, n_atom) - 1]].BreakBond(model[arange[GetHash(nodes, arange, current.element[i], flags, n_atom) - 1]].code, 0);
						model[arange[GetHash(nodes, arange, current.element[i], flags, n_atom) - 1]].BreakBond(model[arange[GetHash(nodes, arange, current.code, flags, n_atom) - 1]].code, 1);	
						model[arange[GetHash(nodes, arange, current.element[i], flags, n_atom) - 1]].BreakBond(model[arange[GetHash(nodes, arange,model[arange[GetHash(node.element[j], flags, n_atom) - 1]].code, 0);
						model[arange[GetHash(nodes, arange,model[arange[GetHash(nodes, arange, current.element[i], flags, n_atom) - 1]].element[j], flags, n_atom) - 1]].BreakBond(model[arange[GetHash(nodes, arange, current.element[i], flags, n_atom) - 1]].code, 1);		
						PrintCompound(model, n_atom);
					}
			}
		}
		for(int i = 0; i < current.element.size(); i++)
		{
			if(!visited[arange[GetHash(nodes, arange, current.element[i], flags, n_atom)]])
			{
				queue.push_back(model[arange[GetHash(nodes, arange, current.element[i], flags, n_atom)]]);
				visited[arange[GetHash(nodes, arange, current.element[i], flags, n_atom)]] = true;
			}
		}
	}
	PrintCompound(model, n_atom);
}

int main()
{
	int n;
	char ch[][3] = {"o1", "c1"};
	
	printf("Enter the number of atoms in your model: ");
	scanf("%d", &n);
	compound model[n];
	int nodes[n], arange[n];
	char atoms[4], element1[4], element2[4];
	for(int i = 0; i < n; i++)
	{
		printf("Enter atom %d : ", i+1);
		cin>>atoms;
		model[i].CreateAtom(atoms);
		nodes[i] = model[i].code;
		arange[i] = i;
	}
	vector<vector<int> > flags;
	flags = SetFlags(nodes, arange, n);
	//PrintArray(nodes, n);
	//PrintArray(arange, n);
	//PrintFlags(flags);
	vector<int> _flag;
	vector<vector<int> >::iterator itr2;
	/*int ch = 0;
	cout<<"Want to have in built standard compounds?(y/n):";
	cin>>ch;
	if(ch)
	{
		printf("1. CARBON DIOXIDE\n2.SULPHUR TRIOXIDE\n3.BENZENE\n4.AMMONIA\n5.WATER\nPlease select one:");
		cin>>ch;
		switch(ch)
		{
			case 1:model
		}
	}*/
	for(itr2 = flags.begin(); itr2 != flags.end(); itr2++)
	{
		_flag = *itr2;
		//PrintVector(_flag);
	}
	int nbonds = 0;
	int count_bonds = 0, temp_hold = 0;
	printf("Enter the number of bonds: ");
	scanf("%d", &nbonds);
	count_bonds = nbonds;
	for(int i = 0; i < nbonds; i++)
	{
		cout<<"Enter the bonding atoms: ";
		cin>>element1>>element2;
		temp_hold = model[arange[GetHash(nodes, arange, encode(element1), flags, n) - 1]].Bond(element2);
		count_bonds -= temp_hold;
		temp_hold = model[arange[GetHash(nodes, arange, encode(element2), flags, n) - 1]].Bond(element1);
	}
	PrintCompound(model, n);
	resonance(model, nodes, arange, flags, count_bonds, n, nbonds);
	return 0;
}
