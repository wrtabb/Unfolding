#ifndef BinningModels_H
#define BinningModels_H

namespace binningmodels{

	// Standard Drell-Yan binning
	const vector<double> _massbinningTrue = {
		15,
		20,
		25,
		30,
		35,
		40,
		45,
		50,
		55,
		60,
		64,
		68,
		72,
		76,
		81,
		86,
		91,
		96,
		101,
		106,
		110,
		115,
		120,
		126,
		133,
		141,
		150,
		160,
		171,
		185,
		200,
		220,
		243,
		273,
		320,
		380,
		440,
		510,
		600,
		700,
		830,
		1000,
		1500,
		3000
	};

	// Standard Drell-Yan binning with each bin split in half
	const vector<double> _massbinningReco0 = {
		15, 
		17.5, 
		20, 
		22.5, 
		25, 
		27.5, 
		30, 
		32.5, 
		35, 
		37.5,
		40,
		42.5,
		45,
		47.5,
		50,
		52.5,
		55,
		57.5,
		60,
		62,
		64,
		66,
		68,
		70,
		72,
		74,
		76,
		78.5,
		81,
		83.5,
		86,
		88.5,
		91,
		93.5,
		96,
		98.5,
		101,
		103.5,
		106,
		108,
		110,
		112.5,
		115,
		117.5,
		120,
		123,
		126,
		129.5,
		133,
		137,
		141,
		145.5,
		150,
		155,
		160,
		165.5,
		171,
		178,
		185,
		192.5,
		200,
		210,
		220,
		231.5,
		243,
		258,
		273,
		296.5,
		320,
		350,
		380,
		410,
		440,
		475,
		510,
		555,
		600,
		650,
		700,
		765,
		830,
		915,
		1000,
		1250,
		1500,
		2250,
		3000
	};

	// Split one bin (last bin)
	const vector<double> _massbinningReco1 = {
		15,
		20,
		25,
		30,
		35,
		40,
		45,
		50,
		55,
		60,
		64,
		68,
		72,
		76,
		81,
		86,
		91,
		96,
		101,
		106,
		110,
		115,
		120,
		126,
		133,
		141,
		150,
		160,
		171,
		185,
		200,
		220,
		243,
		273,
		320,
		380,
		440,
		510,
		600,
		700,
		830,
		1000,
		1500,
		2250,
		3000
	};

	// Split one bin (lowest bin)
	const vector<double> _massbinningReco2 = {
		15,
		17.5,
		20,
		25,
		30,
		35,
		40,
		45,
		50,
		55,
		60,
		64,
		68,
		72,
		76,
		81,
		86,
		91,
		96,
		101,
		106,
		110,
		115,
		120,
		126,
		133,
		141,
		150,
		160,
		171,
		185,
		200,
		220,
		243,
		273,
		320,
		380,
		440,
		510,
		600,
		700,
		830,
		1000,
		1500,
		3000
	};

	// Split lowest bin and highest bin
	const vector<double> _massbinningReco3 = {
		15,
		17.5,
		20,
		25,
		30,
		35,
		40,
		45,
		50,
		55,
		60,
		64,
		68,
		72,
		76,
		81,
		86,
		91,
		96,
		101,
		106,
		110,
		115,
		120,
		126,
		133,
		141,
		150,
		160,
		171,
		185,
		200,
		220,
		243,
		273,
		320,
		380,
		440,
		510,
		600,
		700,
		830,
		1000,
		1500,
		2250,
		3000
	};

}//end namespace

#endif

