#include "format.cuh"

void time_progress(time_t real_start_t, time_t curr_t, double done_part, string proc_name)
{
	time_t real_t = curr_t - real_start_t;
	double left_t = real_t * (1/done_part - 1);
	int b_i = int(left_t);
	time_t b_t = curr_t + b_i;

	cout << proc_name << " " << 100*done_part << " %          \n"
		 << "time used " << real_t/3600 << ":" << (real_t%3600)/60 << ":" << real_t%60 << "          \n"
		 << "time left " << b_i/3600 << ":" << (b_i%3600)/60 << ":" << b_i%60 << "          \n"
		 << "last save: " << string(ctime(&curr_t))
		 << "finish   : " << string(ctime(&b_t))
		 << "\r"      // goto begin & erase
		 << "\033[A"  // up & erase
		 << "\033[A"  // up & erase
		 << "\033[A"  // up & erase
		 << "\033[A"  // up & erase
		 << "\033[A"; // up & erase
	/*
		 * At first sight it seems nothing will be printed because I print and erase all right away.
		 * But I suppose that in some case \smth affect output-buffer but doesn't cause any screen.reprint()
		 * so everything (empty strings from previous cout) is printed in the beginning of next time print when actual chars are printed
		 * So in fact user sees everything erased just before new information to be printed
		 */
}

template<typename T>
//int stp(T str)
void stp(T str)
{
	cerr << str << endl;
	cin.get();
	//int r;
	//cin >> r;
	//return r;
}

template<typename T>
string vectorToStr(vector<T> v, string sp)
{
        string s="";
        for(int i = 0; i < v.size()-1; ++i) s+= (toString(v[i])+sp);
        return s+toString(v[v.size()-1]);
}
template<typename T>
void printVector(ostream &output, vector<T> v, string sp1, string sp2, string sp3)
{
	int sz = v.size();
    output << sp1;
    for(int i = 0; i < sz-1; ++i) output << v[i] << sp2;
    output << v[sz-1] << sp3;
}
template<typename T>
T sumVector(vector<T> v)
{
        T s=0;
        for(int i = 0; i < v.size(); ++i) s+=v[i];
        return s;
}
template <typename T>
string toString(T val)
{
    std::ostringstream oss;
    oss << val;
    return oss.str();
}
template<typename T>
T fromString(const string& s)
{
  std::istringstream iss(s);
  T res;
  iss >> res;
  return res;
}

string toLower(string s)
{
	char d = 'a'-'A';
	for(int i = 0; i<s.size(); ++i) if((s[i]>='A') && (s[i]<='Z')) s[i]+=d;
	return s;
}
string toUpper(string s)
{
	char d = 'A'-'a';
	for(int i = 0; i<s.size(); ++i) if((s[i]>='a') && (s[i]<='z')) s[i]+=d;
	return s;
}
void printX3(ostream &output, double3 r, string sp1, string sp2, string sp3){ output << toStr(r, sp1, sp2, sp3); }
void printX3(ostream &output, int3 r, string sp1, string sp2, string sp3){ output << toStr(r, sp1, sp2, sp3); }

int CHECK(int n)
{
    if(n){
    	ofstream Fout(global_logFname.c_str(), ios::app);
    	//if(!Fout){ return CANT_OPEN_FILE_FOR_WRITING; }
    	if(!Fout){ return n; }

    	if(n == SAY_IT){
    		Fout << error_handl_string;
    	} else {
    		string _s = "error #" + toString(n) + "                                                 \n" +
    				    "message:                                                \n" +
    				    error_handl_string + "\n";
    		Fout << "\n" << _s << "\n";
    		cout << _s;
    	}
    	Fout.close();
    }

    return n;
}

template <typename T>
void SAY_LOG(T s)
{
	error_handl_string = toString(s);
	CHECK(SAY_IT);
}

string toStr(double3 d, string sp1, string sp2, string sp3){ return sp1 + toString(d.x) + sp2 + toString(d.y) + sp2 + toString(d.z) + sp3; }
string toStr(int3 d, string sp1, string sp2, string sp3){ return sp1 + toString(d.x) + sp2 + toString(d.y) + sp2 + toString(d.z) + sp3; }

vector<double> d3ToV(double3 v){
	double vp[3] = {v.x, v.y, v.z};
	vector<double> vv;
	vv.assign(vp, vp+3);
	return vv;
}
