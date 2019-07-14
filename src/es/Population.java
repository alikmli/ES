package es;

import java.util.Arrays;
import java.util.Comparator;

public class Population {
	Individual[] population;
	private double minVal;
	private double maxVal;
	private int dimen;

	public Population(int popSize) {
		this.population = new Individual[popSize];
	}

	public Population(int popSize, int dimen, double minVal, double maxVal, MutationType type) {
		population = new Individual[popSize];
		this.minVal = minVal;
		this.maxVal = maxVal;
		this.dimen = dimen;

		for (int i = 0; i < popSize; i++) {
			population[i] = new Individual(dimen, minVal, maxVal, type);
		}
	}

	public int size() {
		return population.length;
	}

	public Individual[] getPopulation() {
		return population;
	}

	public void setPopulation(Individual[] population) {
		this.population = population;
	}

	public double getMinVal() {
		return minVal;
	}

	public double getMaxVal() {
		return maxVal;
	}

	public int getDimen() {
		return dimen;
	}

	public int populationSize() {
		return population.length;
	}

	public void setIndividual(int index, Individual item) {
		population[index] = item;
	}

	public Individual getIndividual(int index) {
		return population[index];
	}


	public void sortPop() {

		// sort population based on fitness  
		Arrays.sort(population, new Comparator<Individual>() {

			@Override
			public int compare(Individual o1, Individual o2) {
				long item1 = getRealPart(Double.toString(o1.getFitness()));
				long item2 = getRealPart(Double.toString(o2.getFitness()));
				if (Double.compare(item1, item2) != 0) {
					return item1 > item2 ? 1 : item1 < item2 ? -1 : 0;
				} else {
					return Double.compare(o1.getFitness(), o2.getFitness()) > 0 ? 1
							: Double.compare(o1.getFitness(), o2.getFitness()) < 0 ? -1 : 0;
				}

			}
		});

		// sort the objective part of each individual
		for (Individual items : population) {
			Arrays.sort(items.getVars());
		}


	}

	private int getRealPart(String value) {
		if (value.contains(".")) {
			int index = value.indexOf('.');
			return Integer.valueOf(value.substring(0, index));
		} else {
			return Integer.valueOf(value);
		}
	}

}
