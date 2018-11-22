#include "PolynomialFit.h"

PolynomialFit::PolynomialFit(int order)
{
	setOrder(order);
}

PolynomialFit::~PolynomialFit()
{}

void PolynomialFit::setOrder(int n)
{
	if (n < 0)
		_order = 2;
	else
		_order = n;
}

int PolynomialFit::getOrder(void) const
{
	return _order;
}

void PolynomialFit::getCoefs(TVectorD &pars) const
{
	pars.ResizeTo(_last_coefs);
	pars = _last_coefs;
}

Int_t PolynomialFit::operator ()(std::vector<double> &xs_in, std::vector<double> &ys_in,
	TVectorD &pars_out, double in_x0){
	return (*this)(xs_in, ys_in, 0, xs_in.size(), pars_out, in_x0);
}

Int_t PolynomialFit::operator ()(std::vector<double> &xs_in, std::vector<double> &ys_in,
	int offset, int N_points, TVectorD &pars_out, double in_x0) //only for a part of a vector
{
	if (xs_in.size() != ys_in.size()) {
		std::cout<<"PolynomialFit::operator(): Error: x-y data size mismatch"<<std::endl;
		return -1;
	}
	if ((xs_in.size()-offset) < N_points) {
		std::cout<<"PolynomialFit::operator(): Error: N points is out of range:"<<std::endl;
		std::cout<<"\tx.size()="<<xs_in.size()<<" offset="<<offset<<" N_points="<<N_points<<std::endl;
		return -2;
	}
	if (offset < 0) {
		std::cout<<"PolynomialFit::operator(): Error: offset is out of range:"<<std::endl;
		std::cout<<"\tx.size()="<<xs_in.size()<<" offset="<<offset<<" N_points="<<N_points<<std::endl;
		return -3;
	}
	if (N_points < (_order + 1)) {
		std::cout<<"PolynomialFit::operator(): Error: no enough N points for fit:"<<std::endl;
		std::cout<<"\torder="<<_order<<" N_points="<<N_points<<std::endl;
		return -4;
	}

	TMatrixD mat(N_points, _order + 1);
	for (int col = 0; col < mat.GetNcols(); col++)
		for (int row = 0; row < mat.GetNrows(); row++)
			mat[row][col] = pow(xs_in[offset+ row] - in_x0, col);
	TVectorD Y(N_points);
	for (int row = 0; row < Y.GetNrows(); row++)
		Y[row] = ys_in[offset + row];
	TMatrixD mT(mat);
	mT.T();
	TMatrixD In(mT*mat);//because normal assignment like mat = mT*mat does not work! First resizing must be done.
	In.SetTol(1e-40);
	_last_coefs.ResizeTo(_order + 1);
	_last_coefs = (In.Invert())*mT*Y;
	pars_out.ResizeTo(_last_coefs);
	pars_out = _last_coefs;
	if (pars_out.GetNrows()!=(_order+1)) {
		return -5;
	}
	return 0;
}

//=========================================================

DataVector::DataVector(Int_t fit_order, Int_t N_used_): x0_used(0), isCached(false), cache_n_from(0), cache_n_to(0),
		fitter(fit_order), use_left(false), use_right(false), is_set_left(false),
		is_set_right (false), left_value(DBL_MAX), right_value(DBL_MAX)
{
	N_used = std::max(1, N_used_);
}

DataVector::DataVector(std::vector < double> &xx, std::vector<double> &yy, Int_t fit_order, Int_t N_used_): DataVector(fit_order, N_used_)
{
	initialize(xx, yy, fit_order, N_used_);
}
DataVector::~DataVector() {}

void DataVector::initialize(std::vector < double> &xx, std::vector<double> &yy,  Int_t fit_order, Int_t N_used_)
{
	N_used = std::max(1, N_used_);
	fitter.setOrder(fit_order);
	isCached = false;
	if (xx.size()!=yy.size()) {
		std::cout<<"DataVector::initialize(): Error: x and y data size mismatch!"<<std::endl;
		return;
	}
	xs = xx;
	ys = yy;
	double mem;
	for (long unsigned int i=0, i_end = std::min(xs.size(),ys.size());i!=i_end;++i) { //sorting ys(xs) by xs
		auto min_index = i;
		double min_val = xs[i];
		for (auto j=i;j!=i_end;++j)
			if (xs[j]<min_val) {
				min_val = xs[j];
				min_index = j;
			}
		if (min_index!=i) {
			mem = xs[min_index];
			xs[min_index] = xs[i];
			xs[i] = mem;
			mem = ys[min_index];
			ys[min_index] = ys[i];
			ys[i] = mem;
		}
	}
}

void DataVector::setOrder(Int_t ord)
{	isCached = false;
	fitter.setOrder(ord); }

Int_t DataVector::getOrder(void)
{	return fitter.getOrder(); }

void DataVector::setNused(Int_t N)
{	isCached = false;
	N_used = std::max(1, N); }

Int_t DataVector::getNused(void)
{	return N_used; }

std::vector<double> DataVector::getCoefs(void)
{
	std::vector<double> out;
	TVectorD coefs;
	fitter.getCoefs(coefs);
	out.resize(coefs.GetNrows());
	for (long unsigned int i=0, i_end=out.size(); i!=i_end; ++i)
		out[i] = coefs[i];
	return out;
}

void DataVector::use_leftmost(bool use)
{
	use_left = use;
}
void DataVector::use_rightmost(bool use)
{
	use_right = use;
}
void DataVector::set_leftmost(double val)
{
	is_set_left = true;
	left_value = val;
}
void DataVector::unset_leftmost(void)
{
	is_set_left = false;
}
void DataVector::set_rightmost(double val)
{
	is_set_right = true;
	right_value = val;
}
void DataVector::unset_rightmost(void)
{
	is_set_right = false;
}
void DataVector::set_out_value(double val)
{
	set_leftmost(val);
	set_rightmost(val);
}
void DataVector::unset_out_value(void)
{
	unset_leftmost();
	unset_rightmost();
}

void DataVector::push (double x, double y) //do not disrupt order
{
	for (auto i=xs.begin(), i_end = xs.end(), j=ys.begin(), j_end = ys.end(); (i!=i_end)&&(j!=j_end);++i, ++j) {
		if (x==*i)
			return; //do not insert points with equal x
		if (*i>x) {
			xs.insert(i,x);
			ys.insert(j,y);
			isCached = false;
			return;
		}
	}
	xs.insert(xs.end(),x);
	ys.insert(ys.end(),y);
	isCached = false;
}

void DataVector::push_back (double x, double y)//does not check that the new array is ordered, but faster
{
	xs.push_back(x);
	ys.push_back(y);
	isCached = false;
}

void DataVector::erase (std::size_t n)
{
	xs.erase(xs.begin()+n);
	ys.erase(ys.begin()+n);
	isCached = false;
}

void DataVector::clear (void)
{
	xs.clear();
	ys.clear();
	isCached = false;
}

std::size_t DataVector::size (void)
{	return std::min(xs.size(), ys.size());}

double DataVector::getX(std::size_t n)
{	return xs[n]; }

double DataVector::getY(std::size_t n)
{	return ys[n]; }

double DataVector::operator()(double point)
{
	if (xs.empty()||ys.empty())
		return DBL_MAX;
	if (point < xs.front()) {
		if (use_left)
			return ys.front(); //TODO: maybe add scaling
		if (is_set_left)
			return left_value;
	}
	if (point > xs.back()) {
		if (use_right)
			return ys.back();
		if (is_set_right)
			return right_value;
	}
	int n_min, n_max;
	get_indices(point, n_min, n_max);
	if (isCached)
		if ((n_min == cache_n_from) && (n_max == cache_n_to))
			return calculate(point); //same polynomial is used
	cache_n_from = n_min;
	cache_n_to = n_max;
	isCached = true;
	TVectorD temp;
	Int_t ret_code = fitter(xs, ys, n_min, n_max-n_min+1 /*==N_used*/, temp, x0_used);
	if (0==ret_code)
		return calculate(point);
	if (-5==ret_code) {//matrix inversion failed
		if (x0_used!=xs[(n_min+n_max)/2])
			return (*this)(point, xs[(n_min+n_max)/2]);
	}
	return DBL_MAX;
}

double DataVector::operator()(double point, double x0)
{
	isCached = isCached&&(x0_used == x0);
	x0_used = x0;
	return (*this)(point);
}

//[n_min, n_max] are used, not [n_min,n_max). N_used==n_max-n_min+1>=order+1
void DataVector::get_indices(double x, int &n_min, int &n_max)
{
	if (x < xs.front())
	{
		n_min = 0;
		n_max = N_used-1; //order is set to be matching the number of points
		return;
	}
	if (x > xs.back())
	{
		n_max = xs.size() - 1;
		n_min = n_max - N_used + 1;
		return;
	}
	int out_ = 0;
	while ((x > xs[out_]) && (x > xs[out_+1])) ++out_;
	n_min = out_ - (N_used-1) / 2; //asymmetrical interpolation in the case of odd order.
	n_max = n_min + (N_used-1);
	if (n_min < 0)
	{
		n_min = 0;
		n_max = n_min + (N_used-1);
	}
	if (n_max >= xs.size())
	{
		n_max = xs.size() - 1;
		n_min = n_max - N_used+1;
	}
}

double DataVector::calculate(double x)
{
	TVectorD coefs;
	fitter.getCoefs(coefs);
	Int_t order = coefs.GetNrows();
	double out_ = 0;
	for (int o_O = 0; o_O < order; ++o_O)
		out_ += coefs[o_O]*pow(x-x0_used, o_O);
	return out_;
}

