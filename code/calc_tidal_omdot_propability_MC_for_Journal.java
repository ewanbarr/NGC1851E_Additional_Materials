/**
 * Java - Project NGC1851E
 * 
 * Find probability for significant tidal omegadot (only based on angles)
 * 
 * MC sampling of the 4-dim. parameter space: 4 angles
 * 
 * Uses exact expressions for highly eccentric binary, obtained from integrating with Mathematica
 * equations (3.69) with perturbation (3.77) from the Poisson & Will book). 
 * Pulsar parameters from DDGR solution (see paper) 
 * 
 */
package calculations;


import java.io.FileWriter;
import java.io.IOException;
import org.spaceroots.mantissa.random.MersenneTwister; //http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/VERSIONS/JAVA/java.html


public class calc_tidal_omdot_propability_MC_for_Journal {

	private static final long NMC = 10000000000L;

	private static final double DEG2RAD  = Math.PI/180.0;
	private static final double YR       = 365.25*86400.0;
	private static final double C        = 2.99792458E10;
	private static final double GMSUN    = 1.3271244E26;
	private static final double AU       = 14959787070000.0;

	private static final double F0       = 178.700750;
	private static final double F1       = -6.173E-15;
	private static final double F2       =  7E-26;

	private static final double PB       =  7.44790 * 86400;

	private static final double OMDOT    =  0.03468 * DEG2RAD/YR;
	private static final double OMDOTerr =  0.00003 * DEG2RAD/YR;

	private static final double EDOT     = -0.06E-14;
	private static final double EDOTerr  =  1.11E-14;

	private static final double XDOT     =  0.4E-12;
	private static final double XDOTerr  =  0.2E-12;


	/**
	 * MAIN
	 * @param args
	 */
	public static void main(String[] args) {

		double m3 =  1.0; // [Msun]
		double v3 = 10.0; // [km/s]

		String path = "MCresults.dat";

		System.out.println("\n m3 (Msun):       \t" + m3);
		System.out.println(  " v3 (km/s):       \t" + v3);

		//Define sufficiently large range for the distance to the third mass.
		//(distance is a free parameter to maximize omdot for given angles and constraints)
		double r3max = 103*Power(m3, 1.0/3.0); // maximum possible omdot = OMDOTerr
		double r3min = r3max/10.0;             // maximum possible omdot = 1000*OMDOTerr ~ OMDOT

		double n3 = Math.sqrt(GMSUN*(3.9+m3)/Math.pow(r3max*AU, 3.0));
		double P3 = 2*Math.PI/n3;

		//Calculate values for circular motion (just for comparison)

		System.out.println("\n Circular motion values:");

		System.out.println(" maximum r3 (AU):\t" + r3max);
		System.out.println(" P (yr):         \t" + P3/YR);
		System.out.println(" v (km/s):       \t" + r3max*AU*n3/1E5);

		n3 = Math.sqrt(GMSUN*(3.9+m3)/Math.pow(r3min*AU, 3.0));
		P3 = 2*Math.PI/n3;

		System.out.println(" minimum r3 (AU):\t" + r3min);
		System.out.println(" P (yr):         \t" + P3/YR);
		System.out.println(" v (km/s):       \t" + r3min*AU*n3/1E5);

		System.out.println("\n Running MC...");

		MersenneTwister mt = new MersenneTwister((new java.util.Date()).getTime());

		int cntr = 0;
		try {
			FileWriter out = new FileWriter(path);

			for (long iMC = 0; iMC < NMC; iMC++) {

				double cosTheta  = -1.0 + 2.0 * mt.nextDouble();
				double Phi       = 2*Math.PI  * mt.nextDouble();

				double cosThetav = -1.0 + 2.0 * mt.nextDouble();
				double Phiv      = 2*Math.PI  * mt.nextDouble();

				double Theta  = Math.acos(cosTheta);
				double Thetav = Math.acos(cosThetav);

				/*
				 *  find maximum allowed omdot for given angles
				 */			

				double omdotMax = -1E99;
				for (int j = 0; j < 201; j++) {
					double r3 = r3min + (r3max - r3min)*j/200.0;

					//*** frequency derivative boundary conditions ***

					//externally induced F1 / observed F1
					double az    = 0.593008*Cos(Theta) * m3/Power(r3,2);
					double F1acc = -F0*az/C;

					if(Abs(F1acc/F1) > 3) continue;

					//externally induced F2 / observed F2
					double azt = 
							((3.964016E-9 - 1.189205E-8*Power(Cos(Theta),2))*Cos(Thetav) - 
									1.189205E-8*Cos(Phi-Phiv)*Cos(Theta)*Sin(Theta)*Sin(Thetav)) 
							* m3*v3/Power(r3,3);
					double F2acc = -F0*azt/C;

					if(Abs(F2acc/F2) > 3) continue;

					//*** orbital boundary conditions ***

					double eta = 4.1578654E-4 * m3 / Power(r3,3);

					//tidally induced change of eccentricity
					double edot =
							-5.551123 * Power(Cos(Theta),2) + 
							12.092764 * Cos(Phi)*Cos(Theta)*Sin(Theta) - 
							 8.674025 * Cos(Theta)*Sin(Phi)*Sin(Theta) + 
							 2.775561 * Power(Sin(Theta),2) + 
							 6.164007 * Cos(2*Phi)*Power(Sin(Theta),2) + 
							 9.447903 * Cos(Phi)*Sin(Phi)*Power(Sin(Theta),2);
					edot *= (eta/PB);

					if(edot > EDOT+3*EDOTerr || edot < EDOT-3*EDOTerr) continue;

					//tidally induced change of orbital inclination
					double idot = 
							6.155001 * Power(Cos(Theta),2) + 
							7.687440 * Cos(Phi)*Cos(Theta)*Sin(Theta) - 
							3.069228 * Cos(Theta)*Sin(Phi)*Sin(Theta) - 
							9.839474 * Cos(Phi)*Sin(Phi)*Power(Sin(Theta),2) - 
							6.155001 * Power(Sin(Phi),2)*Power(Sin(Theta),2);
					idot *= (eta/PB);

					double xdot = 21.73458*idot;

					if(xdot > XDOT+3*XDOTerr || xdot < XDOT-3*XDOTerr) continue;

					//*** tidally induced change in periastron precession consistent with boundary conditions ***

					double omdot = 
							-6.656531 - 
							 0.057323 * Power(Cos(Theta),2) + 
							13.799088 * Cos(Phi)*Cos(Theta)*Sin(Theta) + 
							26.681274 * Cos(Theta)*Sin(Phi)*Sin(Theta) + 
							10.013457 * Power(Sin(Theta),2) - 
							10.867229 * Cos(2*Phi)*Power(Sin(Theta),2) + 
							23.359612 * Cos(Phi)*Sin(Phi)*Power(Sin(Theta),2);					
					omdot *= (eta/PB);

					if(omdotMax < Abs(omdot)) 
						omdotMax = Abs(omdot);


				} //end of r12 loop

				/*
				 * Check if the maximum induced omdot is significant, 
				 * i.e. > 2 times the error of observed omdot
				 */

				if(omdotMax < 2*OMDOTerr) continue;

				cntr += 1;

				out.write(Phi + "\t" + cosTheta + "\t" + Phiv + "\t" + cosThetav + "\t" + (omdotMax/OMDOT) + "\n");

			} //end of MC loop

			out.close();
		} catch (IOException exc) {
			System.out.println("ERROR: File " + exc.toString());
		}

		System.out.println("\n NMC:                        " + NMC);
		System.out.println(  " Cntr:                       " + cntr);
		System.out.println(  " Cntr / NMC:                 " + ((double)cntr/NMC));
	}

	
	/*
	 * Mathematica CForm functions
	 */
	private static double Abs(double x) { return Math.abs(x); }
	private static double Cos(double x) { return Math.cos(x); }
	private static double Sin(double x) { return Math.sin(x); }
	private static double Power(double x, double a) { return Math.pow(x,a); }


}
