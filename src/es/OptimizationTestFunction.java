package es;

public class OptimizationTestFunction {

	/**
	 * Recommended variable values are: a = 20, b = 0.2 and c = 2π.
	 * 
	 * @param a
	 * @param b
	 * @param c
	 * @param vars
	 * @return
	 */
	public static double Ackley(double[] vars) {
		double term1 = 0;
		double term2 = 0;
		double a = 20, b = 0.2, c = 2 * Math.PI;
		for (int i = 0; i < vars.length; i++) {
			term1 += Math.pow(vars[i], 2);
			term2 += Math.cos(c * vars[i]);
		}

		double d = 1 / vars.length;

		return -a * Math.exp(-b * Math.sqrt(d * term1)) - Math.exp(d * term2) + a + Math.exp(1);
	}

	/**
	 * Input Domain: The function is usually evaluated on the hypercube xi ∈ [-5.12,
	 * 5.12], for all i = 1, …, d. Global Minimum:0
	 * 
	 * @param vars
	 * @return
	 */

	public static double Rastrigin(double[] vars) {
		double term1 = 0;
		for (int i = 0; i < vars.length; i++) {
			term1 += Math.pow(vars[i], 2) - 10 * Math.cos(2 * Math.PI * vars[i]);
		}

		return 10 * vars.length + term1;
	}

	/**
	 * Input Domain: The function is usually evaluated on the hypercube xi ∈ [-5.12,
	 * 5.12], for all i = 1, …, d. Global Minimum:0
	 * 
	 * @param vars
	 * @return
	 */

	public static double Sphere(double[] vars) {
		double term = 0;

		for (int i = 0; i < vars.length; i++) {
			term += Math.pow(vars[i], 2);
		}

		return term;
	}

	/**
	 * Input Domain: The function is usually evaluated on the hypercube xi ∈ [-5,
	 * 10], for all i = 1, …, d, although it may be restricted to the hypercube xi ∈
	 * [-2.048, 2.048], for all i = 1, …, d. Global Minimum:0
	 * 
	 * @param vars
	 * @return
	 */
	public static double Rosenbrock(double[] vars) {
		double term = 0;

		for (int i = 0; i < vars.length - 1; i++) {
			term += 100 * Math.pow((vars[i + 1] - Math.pow(vars[i], 2)), 2) + Math.pow((vars[i] - 1), 2);
		}

		return term;
	}

	/**
	 * Input Domain: The function is usually evaluated on the hypercube xi ∈ [-500,
	 * 500], for all i = 1, …, d. Global Minimum: 0
	 * 
	 * @param vars
	 * @return
	 */

	public static double Schwefel(double[] vars) {
		double term = 0;
		for (int i = 0; i < vars.length; i++) {
			term += vars[i] * Math.sin(Math.sqrt(Math.abs(vars[i])));
		}

		return (418.9829 * vars.length) - term;
	}
}
