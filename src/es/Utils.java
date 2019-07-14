package es;

import java.security.SecureRandom;

public class Utils {
	private static SecureRandom random = new SecureRandom();

	/**
	 * {@link https://github.com/jlizier/jidt/blob/master/java/source/infodynamics/utils/commonsmath3/distribution/UniformIntegerDistribution.java}
	 * 
	 * @param upper
	 * @param lower
	 * @return
	 */
	public static double unformIntRandom(double upper, double lower) {
		final int max =(int) (upper - lower) + 1;
		if (max <= 0) {
			while (true) {
				final int r = random.nextInt();
				if (r >= lower && r <= upper) {
					return r;
				}
			}
		} else {
			return lower +(upper - lower)* random.nextDouble();
		}
	}
	
	/**
	 * 
	 * @return get a random number in a normal distribution
	 */

	public static double unformrealRandom() {

		return random.nextGaussian();
	}
	
	public static double getRandomValue() {
		return random.nextDouble();
	}

}
