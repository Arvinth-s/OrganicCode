//OM NAMO NARAYANA
#include<iostream>
#include<stdio.h>
#include<vector>
#include<list>
#include<string.h>


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

//This will return the index of the least number that is greater than or equal to that number and whether that series is continuous and the number of elements in that sequence
vector<int> Seq(int *arr, int ele, int start)
{
	//cout<<"Entered Seq\n";
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
	//cout<<"Leaving Seq\n";
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
	//cout<<"Entered SetFlags\n";
	int start = 0;
	PrintArray(arr, n);
	QuickSort(arr, arange, start, n-1);
	PrintArray(arr, n);
	vector<int*> flag;
	vector<int> temp_vector;
	vector<vector<int> > index;
	for(int i = 0; i < n; i++)
	{
		index.push_back(Seq(arr, arr[i], i));
		temp_vector = index.back();
		i += temp_vector[1];
	}
	PrintVector(index.back());
	//cout<<"outside for\n";
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

void PrintFlags(vector<int*> flags)
{
	cout<<"Printing flags\n";
	vector<int*>::iterator itr;
	int *temp;
	for(itr = flags.begin(); itr!= flags.end(); itr++)
	{
		temp = *itr;
		printf("%d\t%d\t%d\n", temp[0], temp[1], temp[2]);
	}
}

void swap(int *arr, int i, int j)
{
	//cout<<"Entered swap\n";
	int temp = arr[i];
	arr[i] = arr[j];
	arr[j] = temp;
}

//end is an included index
void QuickSort(int *arr, int *arange, int start, int end)
{
	//cout<<"Entered QuickSort\n";
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
	//cout<<"Entered Search\n";
	int sum = 0;
	vector<vector<int> >::iterator itr;
	vector<int> temp_vector;
	//cout<<"Goin to for loop\n";
	for(itr = flags.begin(); itr!=flags.end();itr++)
	{
		//cout<<"inside for loop\n";
		temp_vector = *itr;
		if(temp_vector[2]/100 > ele/100)
		{
			itr--;
			temp_vector = *itr;
			sum -= temp_vector[1]+1;
			if(!temp_vector[0])
			{
				//cout<<"leaving search with value<<sum<<"<<sum+ele%100 - 1<<"\tsum"<<sum<<"\n";
				return (sum + ele%100);
			}
			else
			{
				for(int k = 0; k< temp_vector[1]+1; k++)
				{
					if(arr[sum+k] == ele)
					{
						//cout<<"leaving search with sum"<<sum+k<<"\n";
						return sum+k;
					}
				}
			}
		}
		else if(temp_vector[2]/100 == ele/100)
		{
			if(!temp_vector[0])
			{
				//cout<<"leaving search with value<<sum<<"<<sum+ele%100 - 1<<"\tsum"<<sum<<"\n";
				return (sum + ele%100);
			}
			else
			{
				for(int k = 0; k< temp_vector[1]+1; k++)
				{
					if(arr[sum+k] == ele)
					{
						//cout<<"leaving search with sum"<<sum+k<<"\n";
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
		//cout<<"Entered GetHash\n";
		//search the position of code in the flags
		int temp_hold = Search(arr, code, flags);
		//cout<<"Finished Search\n";
		if(!temp_hold)
		{
			printf("Something went wrong. Please check the input. Terminating the Program\n");
			return -1;
		}
		else
		{
			//cout<<temp_hold<<"working\n";
			return temp_hold;
		}
	}

void Instruction()
{
	char ch;
	printf("Do you want to read only the rules(n) or understand the model(y): ");
	scanf("%c", &ch);
	if(ch != 'n')
	{
		printf("U should enter the number of atoms in your compounds followed by name of each atom\n. The name of the each atom should be of unique\n");
		printf("The name of each atom should be of the form XXY where XX represents its chemical formula and y its index. Note the index of same element should be different");
		printf("\nFor example if we have 2 Carbons in out model we will represent one as C1 and other as C2\n");
		printf("Then when 'type the number of bonds' msg prompts do the same and add each bond by entering the bond between each compound using the convention given above\n");
		printf("U r ready to go\n\n");		
		return;
	}
	
	else
	{
		printf("READY?!\n\n");
		return;
	}
}


struct keys{
	char _keys[];
	int num;
};
//This function is supposed to return the index of the carbon

//Gives a unique code for every atom, but we cann't use this code directly. It has values in hundreds/thousands.
//So we need a hashtable to store these compounds for quick retrivel

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
		
	char halogens[][3] = {"F", "Cl", "Br", "I"};
	char NFam[][3] = {"N", "P", "As", "Sb", "Bi"};
	char OFam[][3] = {"O", "S", "Se", "Te", "Po"};
	char BFam[][3] = {"B", "Al"};
	char carbon[][3] = {"C"};
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
		return 4000+100+pos;
	
}

//This class contains information of each atom
//Each atom in the compound has its object in compound class, i.e., no of atoms = no. of objects
//V represents the valency of the compounds, since each object represents an atom we can set the maximum no of bonds for that atom

struct compound{
	static int count;
	int key;//index of the object of the structure corresponding to this atom
	char atom[4];//Chemical formula for that atom
	int V;//Maximum possible valency of that atom
	int state;//Oxidation state
	int code;//integer identity
	list<int> element;//Contains the list of bonded atoms' code
	int electronegativity;//Electronegativity of that atom
	int filled;
	void Bond(char batom[]);//To form a bond
	void BondBreak(char batom[]);//To break a bond
	void CreateAtom(char atom[]);//constructor
	int CalcElectroNegativity(int code);
	//int CalcElectroNegativity();
};

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
	(V > 5)?V = 8 - V :V = V;
	this->state = 0;
	this->electronegativity = CalcElectroNegativity(0);
}


void PrintCompound(struct compound *com, int n)
{
	printf("Printing the compound\n");
	for(int i = 0; i<n; i++)
	{
		cout<<"state: "<<com[i].state<<"\tatom:"<<com[i].atom<<"\telectronegativity:"<<com[i].electronegativity<<"\t valency:"<<com[i].V<<endl;
	}
}


//This adds the bonds between the atoms and updates the electronegativity of the atom
void compound::Bond(char batom[])
{
	if(this->element.size() >= V)
	{
		printf("The compound has reached its maximum valency\n");
		return;
	}
	cout<<"electronegativity: "<<CalcElectroNegativity(encode(batom))<<"\t\telectronegativity of same: "<<this->electronegativity<<endl;
	this->element.push_back(encode(batom));
	if (this->electronegativity > CalcElectroNegativity(encode(batom)))
		this->state = this->state-1;
	else if(this->electronegativity < CalcElectroNegativity(encode(batom)))
		this->state = this->state + 1;
	printf("Both have same electronegativity\n");
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
		cin>>atoms;
		model[i].CreateAtom(atoms);
		nodes[i] = model[i].code;
		arange[i] = i;
	}
	PrintArray(nodes, n);
	vector<vector<int> > flags;
	flags = SetFlags(nodes, arange, n);
	vector<int> _flag;
	vector<vector<int> >::iterator itr2;
	for(itr2 = flags.begin(); itr2 != flags.end(); itr2++)
	{
		_flag = *itr2;
		PrintVector(_flag);
	}
	int nbonds = 0;
	printf("Enter the number of bonds\n");
	scanf("%d", &nbonds);
	for(int i = 0; i < nbonds; i++)
	{
		
		cin>>element1>>element2;
		model[GetHash(nodes, arange, encode(element1), flags, n) - 1].Bond(element2);
		cout<<"model1 working\n";
		model[GetHash(nodes, arange, encode(element2), flags, n) - 1].Bond(element1);
		cout<<"model2 working\n";
	}
	PrintCompound(model, n);
	cout<<"Working at last\n";
	return 0;
}
