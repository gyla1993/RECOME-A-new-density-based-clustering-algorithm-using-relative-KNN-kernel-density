#include "head.h"

using namespace std;

struct ELE {   
	int id, pre;
	double w;
	ELE() {}
	ELE(int i_, int p_, double w_) : id(i_), pre(p_), w(w_) {};
	friend bool operator < (const ELE &a, const ELE &b) {
		return a.w < b.w;
	}
};

//Union-Find set 
struct FS {  
	int *fa, *cnt;
	void init(int n) {
		fa = new int[n];
		cnt = new int[n];
		for (int i = 0; i < n; ++i) {
			fa[i] = i;
			cnt[i] = 1;
		}
	}
	int find(int u) {
		return u == fa[u] ? u : (fa[u] = find(fa[u]));
	}
	bool Union(int u, int v) {
		u = find(u);
		v = find(v);
		if (u == v) return false;
		fa[v] = u;
		bool ret = (cnt[u] > 0 && cnt[v] > 0);
		cnt[u] += cnt[v];
		return ret;
	}
	void destroy() {
		delete[] fa;
		delete[] cnt;
	}
};

//FJDD algorithm, which can extract all jump discontinuity values of parameter alpha, stored in alphaArr
void getAlpha(int n, G &sg, vector<double> &alphaArr, double *rhoStar) {
	bool *vis = new bool[n];
	double *val = new double[n];
	vector<ELE> eArr(n - 1);
	eArr.clear();
	priority_queue<ELE> que;

	for (int j = 0; j < n; ++j) {
		val[j] = 0;
		vis[j] = false;
	}
	for (int i = 0; i < n; ++i) {
		if (vis[i] || rhoStar[i] < 1 - 1e-8) continue;
		que.push(ELE(i, -1, rhoStar[i]));
		val[i] = rhoStar[i];
		while (!que.empty()) {
			ELE t = que.top(), t1;
			que.pop();
			if (vis[t.id]) continue;
			val[t.id] = true;
			if (t.pre >= 0) {
				eArr.push_back(t);
			}
			for (E *j = sg.h[t.id]; j; j = j->next) {
				if (vis[j->v]) continue;
				double tw = min(rhoStar[t.id], rhoStar[j->v]);
				if (tw <= val[j->v]) continue;
				val[j->v] = tw;
				que.push(ELE(j->v, t.id, tw));
			}
		}
	}
	sort(eArr.begin(), eArr.end());
	reverse(eArr.begin(), eArr.end());
	FS fs;
	fs.init(n);
	for (int i = 0; i < n; ++i)
		if (rhoStar[i] < 1)
			fs.cnt[i] = 0;
	alphaArr.clear();
	for (vector<ELE>::iterator it = eArr.begin(); it != eArr.end(); ++it) {
		if (fs.Union(it->pre, it->id)) alphaArr.push_back(it->w);
	}
	alphaArr.push_back(0);
	reverse(alphaArr.begin(), alphaArr.end());
	double tt = 100;
	for (vector<double>::iterator it = alphaArr.begin(); it != alphaArr.end(); ++it)
		*it = ceil(*it*tt) / tt;

	//erase duplicate values 
	alphaArr.erase(unique(alphaArr.begin(), alphaArr.end()), alphaArr.end());
	fs.destroy();
	delete[] val;
	delete[] vis;
}