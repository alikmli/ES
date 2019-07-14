package es;

import java.util.Scanner;

public class ESOpt {

	private final static int RHO_SIZE = 2;
	private final static int MU_SIZE = 8;
	private final static int LANDA_SIZE = 56;
	private final static int DIMEN = 5;
	private static double minValue = 0, maxValue = 0;
	private static final int MAX_GEN = 1000;

	public static void main(String[] args) {

		System.out.println("WelCome ....");
		System.out.println("Witch Function you wish to utilized !");
		int counter = 1;
		for (ObjectiveFunctions item : ObjectiveFunctions.values())
			System.out.println((counter++) + " : " + item.name());

		Scanner input = new Scanner(System.in);
		int index = input.nextInt();

		ObjectiveFunctions eval_func = setBorder(index);
		int generation = 1;

		Population pop = new Population(MU_SIZE, DIMEN, minValue, maxValue, MutationType.UNCORRELATED_TYPE1);
		EvaluationSt ES = new EvaluationSt(MU_SIZE, LANDA_SIZE, RHO_SIZE, minValue, maxValue, DIMEN, pop,
				MutationType.UNCORRELATED_TYPE1, eval_func);
		while (generation <= MAX_GEN) {
			ES.recombination(RecombinationType.DOMINANT, RecombinationType.INTERMEDIATE);
			ES.mutation();
			ES.Selection(SelectionType.COMMA,generation);

			System.out.println("\n\n");

			for (Individual item : pop.getPopulation())
				System.out.println(item + "   F = " + item.getFitness());

			generation++;
			System.out.println(ES.getGbest());
		}
	
		System.out.println();
		input.close();

		
	}

	private static ObjectiveFunctions setBorder(int index) {
		switch (index) {
		case 1:
			minValue = -32;
			maxValue = 32;
			return ObjectiveFunctions.ACKLEY;
		case 2:
			minValue = -5.12;
			maxValue = 5.12;
			return ObjectiveFunctions.RASTRIGIN;
		case 3:
			minValue = -5.12;
			maxValue = 5.12;
			return ObjectiveFunctions.SPHERE;
		case 4:
			minValue = -5;
			maxValue = 10;
			return ObjectiveFunctions.ROSENBROCK;
		case 5:
			minValue = -500;
			maxValue = 500;
			return ObjectiveFunctions.SCHWEFE;
		}

		return null;
	}
}
