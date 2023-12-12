#ifndef LS_H_
#define LS_H_
#include <iostream>
#include <limits>
#include <stdint.h>
template <typename T>
class LarssonSadakane
{
	T *I,
		*V,
		r,
		h;
#define LS_KEY(p) (V[*(p) + (h)])
#define LS_SWAP(p, q) (tmp = *(p), *(p) = *(q), *(q) = tmp)
#define LS_SMED3(a, b, c) (LS_KEY(a) < LS_KEY(b) ? (LS_KEY(b) < LS_KEY(c) ? (b) : LS_KEY(a) < LS_KEY(c) ? (c)  \
																										: (a)) \
												 : (LS_KEY(b) > LS_KEY(c) ? (b) : LS_KEY(a) > LS_KEY(c) ? (c)  \
																										: (a)))
	inline void update_group(T *pl, T *pm)
	{
		T g;
		g = (T)(pm - I);
		V[*pl] = g;
		if (pl == pm)
			*pl = -1;
		else
			do
				V[*++pl] = g;
			while (pl < pm);
	}
	inline void select_sort_split(T *p, T n)
	{
		T *pa, *pb, *pi, *pn;
		T f, v, tmp;
		pa = p;
		pn = p + n - 1;
		while (pa < pn)
		{
			for (pi = pb = pa + 1, f = LS_KEY(pa); pi <= pn; ++pi)
				if ((v = LS_KEY(pi)) < f)
				{
					f = v;
					LS_SWAP(pi, pa);
					pb = pa + 1;
				}
				else if (v == f)
				{
					LS_SWAP(pi, pb);
					++pb;
				}
			update_group(pa, pb - 1);
			pa = pb;
		}
		if (pa == pn)
		{
			V[*pa] = (T)(pa - I);
			*pa = -1;
		}
	}
	inline T choose_pivot(T *p, T n)
	{
		T *pl, *pm, *pn;
		T s;
		pm = p + (n >> 1);
		if (n > 7)
		{
			pl = p;
			pn = p + n - 1;
			if (n > 40)
			{
				s = n >> 3;
				pl = LS_SMED3(pl, pl + s, pl + s + s);
				pm = LS_SMED3(pm - s, pm, pm + s);
				pn = LS_SMED3(pn - s - s, pn - s, pn);
			}
			pm = LS_SMED3(pl, pm, pn);
		}
		return LS_KEY(pm);
	}
	inline void sort_split(T *p, T n)
	{
		T *pa, *pb, *pc, *pd, *pl, *pm, *pn;
		T f, v, s, t, tmp;
		if (n < 7)
		{
			select_sort_split(p, n);
			return;
		}
		v = choose_pivot(p, n);
		pa = pb = p;
		pc = pd = p + n - 1;
		while (1)
		{
			while (pb <= pc && (f = LS_KEY(pb)) <= v)
			{
				if (f == v)
				{
					LS_SWAP(pa, pb);
					++pa;
				}
				++pb;
			}
			while (pc >= pb && (f = LS_KEY(pc)) >= v)
			{
				if (f == v)
				{
					LS_SWAP(pc, pd);
					--pd;
				}
				--pc;
			}
			if (pb > pc)
				break;
			LS_SWAP(pb, pc);
			++pb;
			--pc;
		}
		pn = p + n;
		if ((s = (T)(pa - p)) > (t = (T)(pb - pa)))
			s = t;
		for (pl = p, pm = pb - s; s; --s, ++pl, ++pm)
			LS_SWAP(pl, pm);
		if ((s = (T)(pd - pc)) > (t = (T)(pn - pd - 1)))
			s = t;
		for (pl = pb, pm = pn - s; s; --s, ++pl, ++pm)
			LS_SWAP(pl, pm);
		s = (T)(pb - pa);
		t = (T)(pd - pc);
		if (s > 0)
			sort_split(p, s);
		update_group(p + s, p + n - t - 1);
		if (t > 0)
			sort_split(p + n - t, t);
	}
	inline void bucketsort(T *x, T *p, T n, T k)
	{
		T *pi, i, c, d, g;
		for (pi = p; pi < p + k; ++pi)
			*pi = -1;
		for (i = 0; i <= n; ++i)
		{
			x[i] = p[c = x[i]];
			p[c] = i;
		}
		for (pi = p + k - 1, i = n; pi >= p; --pi)
		{
			d = x[c = *pi];
			x[c] = g = i;
			if (d == 0 || d > 0)
			{
				p[i--] = c;
				do
				{
					d = x[c = d];
					x[c] = g;
					p[i--] = c;
				} while (d == 0 || d > 0);
			}
			else
				p[i--] = -1;
		}
	}
	inline T transform(T *x, T *p, T n, T k, T l, T q)
	{
		T b, c, d, e, i, j, m, s;
		T *pi, *pj;
		for (s = 0, i = k - l; i; i >>= 1)
			++s;
		e = std::numeric_limits<T>::max() >> s;
		for (b = d = r = 0; r < n && d <= e && (c = d << s | (k - l)) <= q; ++r)
		{
			b = b << s | (x[r] - l + 1);
			d = c;
		}
		m = (((T)1) << (r - 1) * s) - 1;
		x[n] = l - 1;
		if (d <= n)
		{
			for (pi = p; pi <= p + d; ++pi)
				*pi = 0;
			for (pi = x + r, c = b; pi <= x + n; ++pi)
			{
				p[c] = 1;
				c = (c & m) << s | (*pi - l + 1);
			}
			for (i = 1; i < r; ++i)
			{
				p[c] = 1;
				c = (c & m) << s;
			}
			for (pi = p, j = 1; pi <= p + d; ++pi)
				if (*pi)
					*pi = j++;
			for (pi = x, pj = x + r, c = b; pj <= x + n; ++pi, ++pj)
			{
				*pi = p[c];
				c = (c & m) << s | (*pj - l + 1);
			}
			while (pi < x + n)
			{
				*pi++ = p[c];
				c = (c & m) << s;
			}
		}
		else
		{
			for (pi = x, pj = x + r, c = b; pj <= x + n; ++pi, ++pj)
			{
				*pi = c;
				c = (c & m) << s | (*pj - l + 1);
			}
			while (pi < x + n)
			{
				*pi++ = c;
				c = (c & m) << s;
			}
			j = d + 1;
		}
		x[n] = 0;
		return j;
	}
public:
	void suffixsort(T *x, T *p, T n, T k, T l)
	{
		T *pi, *pk;
		T i, j, s, sl;
		V = x;
		I = p;
		if (n >= k - l)
		{
			j = transform(V, I, n, k, l, n);
			bucketsort(V, I, n, j);
		}
		else
		{
			transform(V, I, n, k, l, std::numeric_limits<T>::max());
			for (i = 0; i <= n; ++i)
				I[i] = i;
			h = 0;
			sort_split(I, n + 1);
		}
		h = r;
		while (*I >= -n)
		{
			pi = I;
			sl = 0;
			do
			{
				if ((s = *pi) <= 0 && (s = *pi) != 0)
				{
					pi -= s;
					sl += s;
				}
				else
				{
					if (sl)
					{
						*(pi + sl) = sl;
						sl = 0;
					}
					pk = I + V[s] + 1;
					sort_split(pi, (T)(pk - pi));
					pi = pk;
				}
			} while (pi <= I + n);
			if (sl)
				*(pi + sl) = sl;
			h = 2 * h;
		}
		for (i = 0; i <= n; ++i)
			I[V[i]] = i;
	}
};
#endif
