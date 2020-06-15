#include "head.h"

using namespace std;

struct CMP {
	double *val;
	CMP() {}
	CMP(double *v_) : val(v_) {}
	bool operator()(int a, int b) {
		return val[a] < val[b];
	}
};

struct CMP1 {
	double *val;
	CMP1() {}
	CMP1(double *v_) : val(v_) {}
	bool operator()(int a, int b) {
		return val[a] > val[b];
	}
};

struct PO {  //data point structure 
	double co[MAXD];       //features of the data point
	friend bool operator < (const PO &a, const PO &b) {
		return a.co[0] == b.co[0] ? a.co[1] < b.co[1] : a.co[0] < b.co[0];
	}
} po[MAXN];

double rho[MAXN], rhoStar[MAXN];     // NKD, RNKD
int output[MAXN], groudtruth[MAXN];  // label outputted by RECOME, groundtruth label
vector<double> alphaArr;  //collection of jump discontinuity values of parameter alpha
int ind_KNN[MAXN][MAXK];   //index of KNN of each point

//some auxiliary variables
bool vis[MAXN];
set<int> KNNSet[MAXN];

template<typename T>
inline bool checkmax(T &a, const T &b) {
	return b > a ? ((a = b), true) : false;
}

template<typename T>
inline bool checkmin(T &a, const T &b) {
	return b < a ? ((a = b), true) : false;
}

template<typename T>
inline T ABS(T a) {
	return a < 0 ? -a : a;
}

//compute the square of the distance between a and b 
double dis2(const PO &a, const PO &b, int dim) { 
	double ret = 0;
	for (int i = 0; i < dim; ++i) ret += (a.co[i] - b.co[i])*(a.co[i] - b.co[i]);
	return ret;
}

void dfs(int u, int label, const double minRho, G &g) {
	output[u] = label;
	for (E *i = g.h[u]; i; i = i->next)
		if (output[i->v] == -1 && rhoStar[i->v] > minRho) {
			dfs(i->v, label, minRho, g);
		}
}

void dfs_assign(int u, int label, G &g1) {
	output[u] = label;
	for (E *i = g1.h[u]; i; i = i->next)
		dfs_assign(i->v, label, g1);
}

int unionAtomCluster(int n, double minRho, G &g, G &g1) {
	int tC = 0;
	for (int i = 0; i < n; ++i) {
		output[i] = -1;
	}
	for (int i = 0; i < n; ++i)
		if (rhoStar[i] == 1 && output[i] == -1) {
			output[i] = ++tC;
			dfs(i, tC, minRho, g);
		}

	for (int i = 0; i < n; ++i)
		if (rhoStar[i] == 1) {
			dfs_assign(i, output[i], g1);
		}
	for (int i = 0; i < n; ++i) {
		if (output[i] > 0) continue;
		output[i] = 1;
	}
	return tC;
}

//n: #points, dim: #dimensions, K, alpha
//if alpha < 0, then it will call FJDD to extract all jump discontinuity values of parameter alpha
int RECOME(int n, int dim, int K, double alpha) {
	int C;  // #clusters
	double sigma = 0;

	//some preparations
	for (int i = 0; i < n; ++i) {
		int tn = 0;
		KNNSet[i].clear();
		for (int j = 1; j <= K; ++j) {
			KNNSet[i].insert(ind_KNN[i][j]);
			--tn;
		}
		sigma += sqrt(dis2(po[i], po[ind_KNN[i][K]], dim));
	}
	sigma /= n;
	//calculate density(NKD)
	for (int i = 0; i < n; ++i) {
		rho[i] = 0;
		for (int j = 1; j <= K; ++j)
			rho[i] += pow(EC, -sqrt(dis2(po[i], po[ind_KNN[i][j]], dim)) / sigma);
	}
	//calculate relative density(RNKD)
	for (int i = 0; i < n; ++i) {
		double mxr = rho[i];
		for (int j = 1; j <= K; ++j) {
			checkmax(mxr, rho[ind_KNN[i][j]]);
		}
		rhoStar[i] = (rho[i] / mxr);
	}

	//find Higher Density Nearest-neighbor (HDN)
	G atom_tree;
	int *pre = new int[n];
	atom_tree.init(n, n);
	for (int i = 0; i < n; ++i) {
		if (rhoStar[i] == 1) {
			pre[i] = i;
			continue;
		}
		double md = LIM;
		for (int j = 1; j <= K; ++j) {
			int tj = ind_KNN[i][j];
			if (rho[tj] > rho[i] && checkmin(md, dis2(po[i], po[tj], dim)))
				pre[i] = tj;
		}
		atom_tree.add(pre[i], i);
	}
	delete pre;

	// construct &KNN graph
	G KNNg;
	KNNg.init(n, K*n * 2);
	for (int i = 0; i < n; ++i)
		for (int j = 1; j <= K; ++j) {
			int t = ind_KNN[i][j];
			if (t <= i) continue;
			if (KNNSet[t].count(i)) {
				KNNg.add(i, t);
				KNNg.add(t, i);
			}
		}

	FILE *fpOut;
	char path[256];

	if (alpha < 0) {
		alphaArr.clear();
		getAlpha(n, KNNg, alphaArr, rhoStar);  //extract all jump discontinuity values of parameter alpha, stored in alphaArr
	}
	else {
		C = unionAtomCluster(n, alpha, KNNg, atom_tree);  //Union atom clusters
	}
	KNNg.destroy();
	atom_tree.destroy();
	return C;   //return #clusters
}

//compute index of KNN of each point, stored in ind_KNN
void calKNN(int n, int mxK, int dim) {  
	double *td = new double[n];
	int *index = new int[n];
	for (int i = 0; i < n; ++i) {
		int cnt = 0;
		for (int j = 0; j < n; ++j) {
			td[j] = dis2(po[i], po[j], dim);
			if (i == j) continue;
			index[cnt++] = j;
		}
		make_heap(index, index + cnt, CMP1(td));
		for (int j = 0; j < mxK; ++j) {
			pop_heap(index, index + cnt - j, CMP1(td));
			ind_KNN[i][j + 1] = index[cnt - j - 1];
		}
	}
	delete td;
	delete index;
}

//calcualate F measure
void calF(int n, int *bel, int *num, double &PB, double &RB, double &F) {
	int *arr, *ass, *ass1;
	ass = new int[n + 1];
	ass1 = new int[n + 1];
	PB = 0, RB = 0;
	for (int i = 0; i < n; ++i) ass[i] = ass1[i] = 0;
	for (int i = 0; i < n; ++i) {
		++ass[bel[i]];
		++ass1[num[i]];
		if (bel[i] < 1)
			throw "shit1";
		if (num[i] < 1)
			throw "shit2";
	}
	for (int i = 0; i < n; ++i) {
		double t = 0, t1 = 0;
		for (int j = 0; j < n; ++j)
			if (i != j) {
				if (bel[i] == bel[j]) t += num[i] == num[j] ? 1 : 0;
				if (num[i] == num[j]) t1 += bel[i] == bel[j] ? 1 : 0;
			}
		if (ass[bel[i]] > 1)PB += t / (ass[bel[i]] - 1);
		if (ass1[num[i]] >1)RB += t1 / (ass1[num[i]] - 1);
	}
	PB /= n;
	RB /= n;
	F = PB*RB * 2 / (PB + RB);
	delete[] ass;
	delete[] ass1;
}

//calcualate NMI
void calNMI(int n, int *bel, int *num, double &NMI) {
	int X = 0, Y = 0;
	for (int i = 0; i < n; ++i) {
		checkmax(X, bel[i]);
		checkmax(Y, num[i]);
	}
	double *cnt = new double[X*Y];
	double *cnt1 = new double[X];
	double *cnt2 = new double[Y];
	for (int i = 0; i < X; ++i)
		for (int j = 0; j < Y; ++j) {
			int t = i*Y + j;
			cnt[t] = 0;
		}
	for (int i = 0; i < X; ++i) cnt1[i] = 0;
	for (int j = 0; j < Y; ++j) cnt2[j] = 0;
	for (int i = 0; i < n; ++i) {
		int l1 = bel[i] - 1, l2 = num[i] - 1;
		cnt[l1*Y + l2] += 1;
		cnt1[l1] += 1;
		cnt2[l2] += 1;
	}
	double MI = 0;
	for (int i = 0; i < X; ++i)
		for (int j = 0; j < Y; ++j) {
			int t = i*Y + j;
			if (cnt[t] < (1 - 1e-3)) continue;
			MI += (cnt[t] / n) * log((cnt[t] * n) / (cnt1[i] * cnt2[j]));
		}
	double H1 = 0, H2 = 0;
	for (int i = 0; i < X; ++i) {
		if (cnt1[i] < (1 - 1e-3)) continue;
		H1 -= (cnt1[i] / n)*log(cnt1[i] / n);
	}
	for (int j = 0; j < Y; ++j) {
		if (cnt2[j] < (1 - 1e-3)) continue;
		H2 -= (cnt2[j] / n)*log(cnt2[j] / n);
	}
	NMI = MI / ((H1 + H2) / 2);
	delete[] cnt;
	delete[] cnt1;
	delete[] cnt2;
}

char path[256];

/***************************
Warning: you may need to adjust pre-defined constants MAXN, MAXK and MAXD 
(defined in head.h) for dataset with large size or high dimension.  
***************************/
int main() {
	FILE *fpIn;
	sprintf(path, "data/in34.txt");
	if ((fpIn = fopen(path, "r+")) == NULL) {
		printf("Open %s Failed!\n", path);
		system("pause");
		exit(-1);
	}
	int n = 0;   //#points
	int dim, C;  //#dimension, #cluster
	int temp;
	bool addOne = false;
	fscanf(fpIn, "%d%d", &dim, &C);
	while (~fscanf(fpIn, "%lf", po[n].co)) {
		for (int i = 1; i < dim; ++i)
			fscanf(fpIn, "%lf", po[n].co + i);
		fscanf(fpIn, "%d", groudtruth + n);
		if (groudtruth[n] < 1)
			addOne = true;
		++n;
	}
	fclose(fpIn);
	if (addOne) {   //ensure that the groudtruth label begins from 1 
		for (int i = 0; i < n; ++i) ++groudtruth[i];
	}


	//compute index of KNN of each point, stored in ind_KNN
	int mxK = sqrt(n) + 1;
	calKNN(n, mxK, dim);

	int K = sqrt(n);
	double alpha = 0.95;

	int tC = RECOME(n, dim, K, alpha);
	double PB, RB, F;
	calF(n, output, groudtruth, PB, RB, F);
	double NMI;
	calNMI(n, output, groudtruth, NMI);
	printf("NMI: %.3f F value: %.3f #Clusters:%d\n", NMI, F, tC);
	return 0;
}


