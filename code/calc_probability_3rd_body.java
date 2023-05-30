/**
 * Project NGC1851E
 * 
 * Find probability for significant tidal omegadot (only based on angles)
 * 
 * MC sampling of the parameter space: 4 angles
 * 
 * Uses exact expressions for a highly eccentic binary 
 * (equations from Poisson & Will, integrated in Mathematica)
 * 
 */
package calculations;


import org.spaceroots.mantissa.random.MersenneTwister;


public class calc_probability_3rd_body {

	private static final long NMC = 10000000000L;

	private static final double DEG2RAD  = Math.PI/180.0;
	private static final double RAD2DEG  = 180.0/Math.PI;
	private static final double YR       = 365.25*86400.0;
	private static final double GMSUN    = 1.3271244E26;
	private static final double AU       = 14959787070000.0;

	private static final double PB       = 7.44790 * 86400;
	private static final double OMDOT    = 0.03468 * DEG2RAD/YR;
	private static final double OMDOTerr = 0.00003 * DEG2RAD/YR;
	private static final double EDOTerr  = 0.011E-12;
	private static final double XDOT     = 0.439E-12;
	private static final double XDOTerr  = 0.202E-12;


	/**
	 * MAIN
	 * @param args
	 */
	public static void main(String[] args) {

		System.out.println("\n===  Calculate probability for significant tidal periastron advance  ========");

		double m3 = 1.0;  // [Msun]
		double v3 = 10.0; // [km/s]

		double r3max = 103*Power(m3, 1.0/3.0); // maximum possible omdot = OMDOTerr
		double r3min = r3max/10.0;             // maximum possible omdot = 1000*OMDOTerr

		double n3 = Math.sqrt(GMSUN*(3.9+m3)/Math.pow(r3max*AU, 3.0));
		double P3 = 2*Math.PI/n3;

		n3 = Math.sqrt(GMSUN*(3.9+m3)/Math.pow(r3min*AU, 3.0));
		P3 = 2*Math.PI/n3;

		MersenneTwister mt = new MersenneTwister((new java.util.Date()).getTime());

		double omdotMax=-1E99,r0=0,Phi0=0,Theta0=0,Phiv0=0,Thetav0=0,edot0=0,xdotx0=0,f10=0,f20=0;
		double omdotMaxTotal=-1E99,r00=0,Phi00=0,Theta00=0,Phiv00=0,Thetav00=0,edot00=0,xdotx00=0,f100=0,f200=0;

		int cntr = 0;
		for (long iMC = 0; iMC < NMC; iMC++) {

			double cosTheta  = -1.0 + 2.0 * mt.nextDouble();
			double Phi       = 2*Math.PI  * mt.nextDouble();

			double cosThetav = -1.0 + 2.0 * mt.nextDouble();
			double Phiv      = 2*Math.PI  * mt.nextDouble();

			double Theta  = Math.acos(cosTheta);
			double Thetav = Math.acos(cosThetav);

			omdotMax = -1E99;
			for (int j = 0; j < 26; j++) {
				double r3 = r3min + (r3max - r3min) *j/25.0;

				double dfF1 = (572191.3737203047*m3*Cos(Theta))/Power(r3,2);

				if(Abs(dfF1) > 3) continue;

				double ddfF2 = (-3.375538221596109e8*m3*v3*(Cos(Thetav) + 
						Cos(Phi)*Cos(Theta)*Cos(Thetav)*Sin(Theta) + 
						Cos(Phi)*Cos(Phi - Phiv)*Power(Sin(Theta),2)*Sin(Thetav)))/Power(r3,3);

				if(Abs(ddfF2) > 3) continue;

				//--- orbital perturbations

				double eta = 0.000415786545078928 * m3 / Power(r3,3);

				double edot = 
						-6.326190236800033 * Power(Cos(Theta),2) + 
						12.930963225519722 * Cos(Phi)*Cos(Theta)*Sin(Theta) - 
						8.119565862683377  * Cos(Theta)*Sin(Phi)*Sin(Theta) + 
						3.163095118400039  * Power(Sin(Theta),2) + 
						5.768428950085463  * Cos(2*Phi)*Power(Sin(Theta),2) + 
						8.298344157181887  * Cos(Phi)*Sin(Phi)*Power(Sin(Theta),2);
				edot *= (eta/PB);
				
				if(Abs(edot) > 2*EDOTerr)      continue;

				double idot = 
						5.7632211086570795 * Power(Cos(Theta),2) + 
						6.737246916772856  * Cos(Phi)*Cos(Theta)*Sin(Theta) - 
						5.28208491411928   * Cos(Theta)*Sin(Phi)*Sin(Theta) - 
						10.498370575126906 * Cos(Phi)*Sin(Phi)*Power(Sin(Theta),2) - 
						5.7632211086570795 * Power(Sin(Phi),2)*Power(Sin(Theta),2);
				idot *= (eta/PB);

				double xdotx = 17.8493174204697 * idot;

				if(Abs(xdotx) > XDOT+2*XDOTerr) continue;
				
				double omdot = 
						-6.65589094875108 + 
						4.755943120079483  * Power(Cos(Theta),2) + 
						16.839370143668887 * Cos(Phi)*Cos(Theta)*Sin(Theta) + 
						28.106756088133793 * Cos(Theta)*Sin(Phi)*Sin(Theta) + 
						7.605864863086828  * Power(Sin(Theta),2) - 
						8.472770982857119  * Cos(2*Phi)*Power(Sin(Theta),2) + 
						20.4746143862389   * Cos(Phi)*Sin(Phi)*Power(Sin(Theta),2);
				omdot *= (eta/PB);

				if(omdotMax < Abs(omdot)) {
					omdotMax = Abs(omdot);
					r0      = r3;
					Phi0    = Phi;
					Theta0  = Theta;
					Phiv0   = Phiv;
					Thetav0 = Thetav;
					edot0   = edot;
					xdotx0  = xdotx;
					f10     = dfF1;
					f20     = ddfF2;
				}

			}

			//no significant solution found
			if(omdotMax < 2*OMDOTerr) continue;

			cntr += 1;

			if(omdotMaxTotal < omdotMax) {
				omdotMaxTotal = omdotMax;
				r00      = r0;
				Phi00    = Phi0;
				Theta00  = Theta0;
				Phiv00   = Phiv0;
				Thetav00 = Thetav0;
				edot00   = edot0;
				xdotx00  = xdotx0;
				f100     = f10;
				f200     = f20;
			}

		}

		n3 = Math.sqrt(GMSUN*(3.9+m3)/Math.pow(r00*AU, 3.0));
		P3 = 2*Math.PI/n3;

		System.out.println("\n NMC:                        " + NMC);
		System.out.println(  " Cntr:                       " + cntr);
		System.out.println(  " Cntr / NMC:                 " + ((double)cntr/NMC));
		System.out.println("\n max. tidal omdot (deg/yr):  " + omdotMaxTotal * RAD2DEG*YR);
		System.out.println(  " -- '' --        / OMDOT:    " + omdotMaxTotal / OMDOT);
		System.out.println(  " -- '' --        / OMDOTerr: " + omdotMaxTotal / OMDOTerr);
		System.out.println("\n r0 (AU):\t" + r00);
		System.out.println(  " P2 (yr):\t" + P3/YR);
		System.out.println("\n Phi   (deg):\t" + Phi00 *RAD2DEG + "\tTheta   (deg): " + Theta00 *RAD2DEG);
		System.out.println(  " Phi_v (deg):\t" + Phiv00*RAD2DEG + "\tTheta_v (deg): " + Thetav00*RAD2DEG);
		System.out.println("\n edot:  \t" + edot00); 
		System.out.println(  " xdotx: \t" + xdotx00); 
		System.out.println(  " df/F1: \t" + f100); 
		System.out.println(  " ddf/F2:\t" + f200); 

		System.out.println("\n=== done ====================================================================");
	}

	/*
	 * Mathematica functions
	 */
	private static double Abs(double x) { return Math.abs(x); }
	private static double Cos(double x) { return Math.cos(x); }
	private static double Sin(double x) { return Math.sin(x); }
	private static double Power(double x, double a) { return Math.pow(x,a); }


}
