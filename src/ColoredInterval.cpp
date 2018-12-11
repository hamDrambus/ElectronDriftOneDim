#include "ColoredInterval.h"

class ColoredRange;


ColoredInterval::ColoredInterval(long double left, long double right, long double color):
		left_(std::min(left, right)), right_(std::max(left, right)), color_(std::fabs(color))
{}

ColoredRange::ColoredRange(void) {}
ColoredRange::ColoredRange(const ColoredInterval &inter)
{
	arr_.push_back(inter);
}

long int ColoredRange::NumOfIndices(void)
{
	long int N = 0;
	for (int j=0, end_ = arr_.size(); j!=end_; ++j) {
		N+=std::lround((arr_[j].right_-arr_[j].left_)/arr_[j].color_);
		if (j==(end_-1)) {
			N+=1;
		} else {
			if (arr_[j+1].left_!=arr_[j].right_)
				N+=1; //intervals are not adjacent
		}
	}
	return N;
}

long double ColoredRange::Value (long int index)
{
	long int Nsum = 0, Noff = 0, N = 0;
	int j=0, end_ = arr_.size();
	bool closed = false;
	for (; j!=end_; ++j) {
		N = std::lround((arr_[j].right_-arr_[j].left_)/arr_[j].color_);
		if (j==(end_-1)) {
			N+=1;
			closed= true;
		} else {
			closed = false;
			if (arr_[j+1].left_!=arr_[j].right_) {
				N+=1; //intervals are not adjacent
				closed = true;
			}
		}
		Nsum+=N;
		if (index<Nsum) {
			Noff = Nsum - N;
			break;
		}
	}
	if (j==end_)
		return DBL_MAX;
	return arr_[j].left_+(index-Noff)*(arr_[j].right_-arr_[j].left_)/(N- (closed ?  1 : 0));
}

void ColoredRange::Print(std::ostream & str)
{
	for (int j=0, end_ = arr_.size(); j!=end_; ++j) {
		str<<"["<<arr_[j].left_<<"; "<<arr_[j].color_<<"; "<<arr_[j].right_<<"]"<<std::endl;
	}
}

ColoredRange operator+ (ColoredRange l, const ColoredRange& r)
{
	for (int j=0, end_=r.arr_.size(); j!=end_; ++j) {
		l = l + r.arr_[j];
	}
	return l;
}

ColoredRange operator+ (ColoredRange l, const ColoredInterval& r)
{
	std::vector<ColoredInterval> new_l;
	bool is_inside_r = false;
	for (int j=0, end_=l.arr_.size(); j!=end_; ++j) {
		if (is_inside_r) {
			if (r.right_<l.arr_[j].left_) {
				new_l.push_back(ColoredInterval(l.arr_[j-1].right_, r.right_, r.color_));
				new_l.insert(new_l.end(),l.arr_.begin()+j,l.arr_.end());
				break;
			}
			if (r.right_<l.arr_[j].right_) {
				new_l.push_back(ColoredInterval(l.arr_[j-1].right_, l.arr_[j].left_, r.color_));
				new_l.push_back(ColoredInterval(l.arr_[j].left_, r.right_, std::min(r.color_, l.arr_[j].color_)));
				new_l.push_back(ColoredInterval(r.right_, l.arr_[j].right_, l.arr_[j].color_));
				new_l.insert(new_l.end(),l.arr_.begin()+j+1,l.arr_.end());
				break;
			}
			//!(r.right_<l.arr_[j].right_)
			new_l.push_back(ColoredInterval(l.arr_[j-1].right_, l.arr_[j].left_, r.color_));
			new_l.push_back(ColoredInterval(l.arr_[j].left_, l.arr_[j].right_, std::min(r.color_, l.arr_[j].color_)));
		} else {
			if (r.right_<l.arr_[j].left_) {
				new_l.push_back(ColoredInterval(r.left_, r.right_, r.color_));
				new_l.insert(new_l.end(),l.arr_.begin()+j,l.arr_.end());
				break;
			}
			if ((r.right_<l.arr_[j].right_)&&(r.left_<l.arr_[j].left_)) {
				new_l.push_back(ColoredInterval(r.left_, l.arr_[j].left_, r.color_));
				new_l.push_back(ColoredInterval(l.arr_[j].left_, r.right_, std::min(r.color_, l.arr_[j].color_)));
				new_l.push_back(ColoredInterval(r.right_, l.arr_[j].right_, l.arr_[j].color_));
				new_l.insert(new_l.end(),l.arr_.begin()+j+1,l.arr_.end());
				break;
			}
			if ((r.right_<l.arr_[j].right_)&&!(r.left_<l.arr_[j].left_)) {
				new_l.push_back(ColoredInterval(l.arr_[j].left_, r.left_, l.arr_[j].color_));
				new_l.push_back(ColoredInterval(r.left_, r.right_, std::min(r.color_, l.arr_[j].color_)));
				new_l.push_back(ColoredInterval(r.right_, l.arr_[j].right_, l.arr_[j].color_));
				new_l.insert(new_l.end(),l.arr_.begin()+j+1,l.arr_.end());
				break;
			}
			if (!(r.left_<l.arr_[j].right_)) {
				new_l.push_back(ColoredInterval(l.arr_[j].left_, l.arr_[j].right_, l.arr_[j].color_));
			} else { //!(r.right_<l.arr_[j].right_) && (r.left_<l.arr_[j].right_)
				is_inside_r = true;
				if (r.left_<l.arr_[j].left_) {
					new_l.push_back(ColoredInterval(r.left_, l.arr_[j].left_, r.color_));
					new_l.push_back(ColoredInterval(l.arr_[j].left_, l.arr_[j].right_, std::min(r.color_, l.arr_[j].color_)));
				} else {
					new_l.push_back(ColoredInterval(l.arr_[j].left_, r.left_, l.arr_[j].color_));
					new_l.push_back(ColoredInterval(r.left_, l.arr_[j].right_, std::min(r.color_, l.arr_[j].color_)));
				}
			}
		}
		if ((end_-1) == j) { //special conditions for edges
			if (is_inside_r) {
				new_l.push_back(ColoredInterval(l.arr_[j].right_, r.right_, r.color_));
			} else { //(r.left_>=l.arr_[j].right_)
				new_l.push_back(ColoredInterval(r.left_, r.right_, r.color_));
			}
		}
	}
	l.arr_.clear();
	for (int j =0, end_ = new_l.size(); j!=end_; ++j) {
		if (l.arr_.empty()) {
			if (new_l[j].right_>new_l[j].left_)
				l.arr_.push_back(new_l[j]);
			continue;
		}
		if ((l.arr_.back().color_==new_l[j].color_)&&(l.arr_.back().right_==new_l[j].left_))
			l.arr_.back().right_ = new_l[j].right_; //merging
		else {
			if (new_l[j].right_>new_l[j].left_)
				l.arr_.push_back(new_l[j]);
		}
	}
	return l;
}

ColoredRange operator+ (ColoredInterval l, const ColoredRange& r)
{
	ColoredRange left (l);
	return left+r;
}

ColoredRange operator+ (ColoredInterval l, const ColoredInterval& r)
{
	ColoredRange left (l);
	return left+r;
}
