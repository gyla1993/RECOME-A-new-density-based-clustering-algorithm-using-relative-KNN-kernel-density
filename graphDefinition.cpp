#include "head.h"

bool E::operator < (const E &t) {
	return w < t.w;
}

void G::init(int n, int en) {
	h = new E *[n];
	e = new E[en];
	memset(h, 0, sizeof(h[0])*n);
	re = e;
}

void G::add(int u, int v, double w) {
	re->u = u;
	re->v = v;
	re->w = w;
	re->next = h[u];
	h[u] = re++;
}

void G::destroy() {
	delete[] h;
	delete[] e;
}

