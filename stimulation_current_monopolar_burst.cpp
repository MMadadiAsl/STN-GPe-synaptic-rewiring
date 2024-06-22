
//***************************************************************************
//***************************************************************************

	#include <iostream>
	#include <cmath>
	#include <fstream>

	using namespace std;

//***************************************************************************
//***************************************************************************

int main(){

	const float PI = 4. * atan(1.);

	int p, np, nb;
	float t, tp[5], amp_stim, I_stim, delta_stim, freq_stim, intra_burst, T_OFF;

	p = 4;
	np = 0;
	nb = 0;
	T_OFF = 500.0;
	amp_stim = 200.0;
	freq_stim = 15.0;
	delta_stim = 5.0;
	intra_burst = 1.0 / freq_stim;

	ofstream outfile1("output1.txt");

//***************************************************************************
//***************************************************************************

	for (t = 0.0; t <= 3000.0; t = t + 0.01){
	
		tp[np] = (1000.0 * np * intra_burst) + (nb * (T_OFF + 291.7));
	
		if (tp[np] <= t && t <= (tp[np] + delta_stim)) I_stim = amp_stim;
			else I_stim = 0.0;
		if (t > (tp[np] + delta_stim)) np++;
		if (np > p){
			np = 0;
			nb++;
		}

	outfile1 << t << '\t' << I_stim << endl;

	}

//***************************************************************************
//***************************************************************************

	outfile1.close();

return 0;
}

//***************************************************************************
//***************************************************************************