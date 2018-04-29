#include<iostream>
#include<vector>
#include<utility>
#include<numeric>
#include<stack>
#include<algorithm>
#include<random>
#include<cmath>
using namespace std;

unsigned long long BinaryConvertToInteger(const vector<short>& binary)
{
	unsigned long long digit = 0;
	for (size_t i = 0; i < binary.size(); ++i)
	{
		if (binary[i] == 1)
		{
			digit += static_cast<decltype(digit)>(pow(2, i));
		}
	}
	return digit;
}


vector<short> IntegerConvertToBinary(unsigned long long digit)
{
	vector<short> binary;
	while (true)
	{
		binary.emplace_back(static_cast<short>(digit % 2));
		digit = digit / 2;
		if (digit == 0)
			break;
	}

	return binary;
}
double func(double x1, double x2)
{
	return (-2 * x2*x2*x2 + 6 * x2*x2 + 6 * x2 + 10)*sin(log(x1)*exp(x2));
}
template<typename TPopulation>
void Init(TPopulation& popul, size_t SCountGenes, size_t Batch)
{
	uniform_int_distribution<int>uid(0, 1);
	random_device rd;
	for (size_t i = 0; i < Batch; ++i)
	{
		vector<short> temp;
		for (size_t j = 0; j < SCountGenes; ++j)
		{
			temp.emplace_back(uid(rd));
		}
		popul.emplace_back(temp);
	}
}
template<typename TConstraints, typename TBits>
double Decoder(TConstraints& cnstr,TBits bits)
{
	return cnstr.first + BinaryConvertToInteger(bits)*(cnstr.second - cnstr.first) / (pow(2, bits.size()) - 1);
}


int main()
{
	vector<pair<double, double>> constraints; // ограничения на переменные
	vector<int> count_genes; // количество генов 
	vector<double> valuesFunction;
	vector<vector<short>> population;// популяция
	double eps = 100000;
	size_t countVariables = 2;



	constraints.emplace_back(0.5, 1.1);
	constraints.emplace_back(1.0, 4.6);

	for (const auto& x : constraints)
	{
		double lower = log2((x.second - x.first)*eps + 1);
		double upper = log2((x.second - x.first)*eps) + 1;
		count_genes.emplace_back(static_cast<int>((lower + upper) / 2 + 0.5));
	}

	size_t SharedCountGenes = accumulate(count_genes.cbegin(), count_genes.cend(), 0);

	Init(population, SharedCountGenes, 10);
	cout << "Initial population" << endl;

	for (const auto& x : population)
	{
		vector<double> vals;
		auto iterBegin = x.cbegin();
		for (size_t i = 0; i < 2; ++i)
		{
			vector<short> temp(count_genes[i]);
			copy(iterBegin, iterBegin + count_genes[i], temp.rbegin());
			vals.emplace_back(Decoder(constraints[i], temp));
			iterBegin = x.cbegin() + count_genes[i];
		}
		cout << func(vals[0], vals[1]) << endl;
	}

	size_t iteration = 0;

	while (iteration < 10000)
	{

		// Вычисление значения функции пригодности хромосом
		for (const auto& x : population)
		{
			vector<double> vals;
			auto iterBegin = x.cbegin();
			for (size_t i = 0; i < 2; ++i)
			{
				vector<short> temp(count_genes[i]);
				copy(iterBegin, iterBegin + count_genes[i], temp.rbegin());
				vals.emplace_back(Decoder(constraints[i], temp));
				iterBegin = x.cbegin() + count_genes[i];
			}
			valuesFunction.emplace_back(func(vals[0], vals[1]));
		}


		//СЕЛЕКЦИЯ
		vector<double> Probab;
		double SharedFunctionRez = 0;
		auto ItMinValue = min_element(valuesFunction.cbegin(), valuesFunction.cend());
		for (const auto& x : valuesFunction)
		{
			SharedFunctionRez += x - (*ItMinValue);
		}

		//вероятность отбора
		for (size_t i = 0; i < valuesFunction.size(); ++i)
		{
			Probab.emplace_back((valuesFunction[i] - (*ItMinValue)) / SharedFunctionRez);
		}

		//совокупная вероятность
		for (size_t i = 1; i < Probab.size(); ++i)
		{
			Probab[i] += Probab[i - 1];
		}
		// рулетка
		uniform_real_distribution<double> urd(0, 1);
		random_device rd;
		vector<vector<short>>buff_population;

		for (size_t i = 0; i < population.size(); ++i)
		{
			double temp = urd(rd);
			for (size_t j = 0; j < population.size(); j++)
			{
				if (Probab[j] >= temp)
				{
					buff_population.emplace_back((population[j]));
					break;
				}
			}

		}

		if (!buff_population.empty())
			population = move(buff_population);

		//КРОССИНГОВЕР
		size_t PointCross;
		double ProbabilityCross = 0.5;
		stack<size_t> index;
		for (size_t i = 0; i < population.size(); ++i)
		{
			double temp = urd(rd);
			if (temp < ProbabilityCross)
			{
				index.push(i);
			}
			if (((index.size() % 2) == 0) && (!index.empty()))
			{
				uniform_int_distribution<size_t> uid(0, SharedCountGenes);
				PointCross = uid(rd);
				vector<short> temp1 = population[index.top()];
				index.pop();
				vector<short>temp2 = population[index.top()];
				index.pop();
				for (size_t k = PointCross + 1; k < SharedCountGenes; ++k)
				{
					short tmp = temp1[k];
					temp1[k] = temp2[k];
					temp2[k] = tmp;
				}

			}
		}
		while (!index.empty())
		{
			index.pop();
		}

		//МУТАЦИЯ

		double ProbabMutation = 0.01;
		for (size_t i = 0; i < population.size(); ++i)
		{
			for (size_t j = 0; j < SharedCountGenes; ++j)
			{
				if (urd(rd) < ProbabMutation)
				{
					if (population[i][j] == 1)
						population[i][j] = 0;
					else
						population[i][j] = 1;
				}
			}
		}

		++iteration;

		//Reset
		valuesFunction.clear();

	}
	cout << endl << "Childs" << endl;

	for (const auto& x : population)
	{
		vector<double> vals;
		auto iterBegin = x.cbegin();
		for (size_t i = 0; i < 2; ++i)
		{
			vector<short> temp(count_genes[i]);
			copy(iterBegin, iterBegin + count_genes[i], temp.rbegin());
			vals.emplace_back(Decoder(constraints[i], temp));
			iterBegin = x.cbegin() + count_genes[i];
		}
		cout << func(vals[0], vals[1]) << endl;
	}
}
