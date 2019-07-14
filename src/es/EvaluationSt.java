package es;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 * 
 * @author Alikamali Eval. Strategy Main Class
 *
 */
public class EvaluationSt {
	/**
	 * mutation probability
	 */
	private final double Pm = 0.4;

	/**
	 * recombination probability
	 */
	private final double Pr = 0.8;
	private int muSize;
	private int lamdaSize;
	private int rhoO;
	private double minVal;
	private double maxVal;
	private int dimen;
	/**
	 * lamada and mu Populdation
	 */
	private Population muPop, lamdaPop;
	/**
	 * epsolon constant for limiting sigma
	 */
	private final double epsilon = 0.001;
	private ObjectiveFunctions functionEval;
	private GBest gbest = new GBest();
	MutationType muType;

	/**
	 * constructor set the fitness of initial population and set a population for
	 * lamda (Empty)<br />
	 * <b> base on mutation type u set instruction of your individual</b>
	 * 
	 * @param muSize       Population size(parents)
	 * @param lamdaSize    lamda Population size(children)
	 * @param rhoO         ρ factory for number of parent will be selected randomly
	 * @param minVal       min value accourding according eval. function
	 * @param maxVal       max value accourding according eval. function
	 * @param dimen        dimension of the individuals
	 * @param muPop        mu Populdation
	 * @param mutationType type of mutation will be used (uncorrelated with 1 step
	 *                     size and n step size also correlated mutation)
	 * @param functionEval eval. function that will be used
	 */
	public EvaluationSt(int muSize, int lamdaSize, int rhoO, double minVal, double maxVal, int dimen, Population muPop,
			MutationType mutationType, ObjectiveFunctions functionEval) {
		this.muSize = muSize;
		this.lamdaSize = lamdaSize;
		this.rhoO = rhoO;
		this.minVal = minVal;
		this.maxVal = maxVal;
		this.dimen = dimen;
		this.muPop = muPop;
		this.functionEval = functionEval;
		this.muType = mutationType;

		for (int i = 0; i < muPop.size(); i++) {
			muPop.getIndividual(i).setFitness(getFitness(muPop.getIndividual(i).getVars()));
		}

		this.lamdaPop = new Population(lamdaSize, dimen, minVal, maxVal, mutationType);
	}

	public Population getPopulation() {
		return muPop;
	}

	public GBest getGbest() {
		return gbest;
	}

	/**
	 * 
	 * @return a Population with ρ number of member
	 */
	public Population getSampleRhoMu() {
		Population tmp = new Population(rhoO);

		for (int i = 0; i < rhoO; i++) {
			int index = (int) (Utils.unformIntRandom(muSize, 1) - 1);
			tmp.setIndividual(i, muPop.getIndividual(index));
		}

		return tmp;
	}

	/**
	 * get list of specified column in the individuals (sigmas or angles)
	 * 
	 * @param tmpPop witch population
	 * @param index  which column
	 * @param type   type=1 sigmas , types=3 angles
	 * @return list of items
	 */
	//
	public List<Double> getColumnsAt(Population tmpPop, int index, int type) {
		List<Double> numbs = new LinkedList<>();
		if (type == 1) {
			for (Individual item : tmpPop.population) {
				numbs.add(item.getSigmaAt(index));
			}
		} else if (type == 2) {
			for (Individual item : tmpPop.population) {
				numbs.add(item.getAngleAt(index));
			}
		}
		return numbs;
	}

	/**
	 * u can recombine in two way, Dominant and intermediate for each part
	 * 
	 * @param recTypeO type for objective vars
	 * @param recTypeS type for strategy part
	 */
	public void recombination(RecombinationType recTypeO, RecombinationType recTypeS) {
		// Objective Part

		if (recTypeO == RecombinationType.DOMINANT) {
			for (int i = 0; i < lamdaSize; i++) {
				Population tmp = getSampleRhoMu();

				if (Utils.getRandomValue() < Pr) {
					for (int j = 0; j < dimen; j++) {
						int randIndex = (int) (Utils.unformIntRandom(tmp.size(), 1) - 1);
						lamdaPop.getIndividual(i).setValue(j, tmp.getIndividual(randIndex).getValue(j));
					}
					// Strategy Parameter part
					StrategyParamRecombination(tmp, recTypeS, i);
				} else {
					int randIndex = (int) (Utils.unformIntRandom(tmp.size(), 1) - 1);
					lamdaPop.setIndividual(i, tmp.getIndividual(randIndex));
				}
			}
		} else {
			for (int i = 0; i < lamdaSize; i++) {
				Population tmp = getSampleRhoMu();

				if (Utils.getRandomValue() < Pr) {
					for (int j = 0; j < dimen; j++) {
						lamdaPop.getIndividual(i).setValue(j, avg(tmp.getPopulation(), j));
					}
					// Strategy Parameter part
					StrategyParamRecombination(tmp, recTypeS, i);

				} else {
					int randIndex = (int) (Utils.unformIntRandom(tmp.size(), 1) - 1);
					lamdaPop.setIndividual(i, tmp.getIndividual(randIndex));
				}
			}
		}

	}

	/**
	 * recombine strategy part base on mutation type(your structure is different)
	 * 
	 * @param tmp      random population
	 * @param recTypeS type of recombination for strategy part
	 * @param i        index of lamda population element
	 */
	private void StrategyParamRecombination(Population tmp, RecombinationType recTypeS, int i) {
		if (recTypeS == RecombinationType.DOMINANT) {
			if (muType == MutationType.UNCORRELATED_TYPE1) {
				int randIndex = (int) (Utils.unformIntRandom(tmp.size(), 1) - 1);
				lamdaPop.getIndividual(i).setSigma(tmp.getIndividual(randIndex).getSigma());
			} else if (muType == MutationType.UNCORRELATED_TYPE2) {
				for (int k = 0; k < dimen; k++) {
					int randIndex = (int) (Utils.unformIntRandom(tmp.size(), 1) - 1);
					lamdaPop.getIndividual(i).setSigmaAt(k, tmp.getIndividual(randIndex).getSigmaAt(k));
				}
			} else if (muType == MutationType.CORRELATED) {
				for (int k = 0; k < dimen; k++) {
					int randIndex = (int) (Utils.unformIntRandom(tmp.size(), 1) - 1);
					lamdaPop.getIndividual(i).setSigmaAt(k, tmp.getIndividual(randIndex).getSigmaAt(k));
				}
				double tmpLen = dimen * (dimen - 1) / 2;
				for (int k = 0; k < tmpLen; k++) {
					int randIndex = (int) (Utils.unformIntRandom(tmp.size(), 1) - 1);
					lamdaPop.getIndividual(i).setAngleAt(k, tmp.getIndividual(randIndex).getAngleAt(k));
				}
			}
		} else {
			if (muType == MutationType.UNCORRELATED_TYPE1) {
				double sum = 0, avg = 0;
				for (int j = 0; j < tmp.size(); j++) {
					sum += tmp.getIndividual(j).getSigma();
				}
				avg = sum / tmp.size();
				lamdaPop.getIndividual(i).setSigma(avg);
			} else if (muType == MutationType.UNCORRELATED_TYPE2) {
				for (int k = 0; k < dimen; k++) {
					List<Double> tmpValues = getColumnsAt(tmp, k, 1);
					double sum = 0, avg = 0;
					for (Double items : tmpValues) {
						sum += items;
					}
					avg = sum / tmpValues.size();
					lamdaPop.getIndividual(i).setSigmaAt(k, avg);
				}
			} else if (muType == MutationType.CORRELATED) {
				// rec. angle
				for (int k = 0; k < dimen; k++) {
					List<Double> tmpValues = getColumnsAt(tmp, k, 1);
					double sum = 0, avg = 0;
					for (Double items : tmpValues) {
						sum += items;
					}
					avg = sum / tmpValues.size();
					lamdaPop.getIndividual(i).setSigmaAt(k, avg);
				}
				double tmpLen = dimen * (dimen - 1) / 2;
				for (int k = 0; k < tmpLen; k++) {
					List<Double> tmpValues = getColumnsAt(tmp, k, 2);
					double sum = 0, avg = 0;
					for (Double items : tmpValues) {
						sum += items;
					}
					avg = sum / tmpValues.size();
					lamdaPop.getIndividual(i).setAngleAt(k, avg);
				}
			}
		}
	}

	/**
	 * 
	 * @param inds target individuals
	 * @param col  which column from it
	 * @return average of specified column from object part all individuals
	 */
	public double avg(Individual[] inds, int col) {
		double avg = 0, sum = 0;

		for (int i = 0; i < rhoO; i++) {
			sum += inds[i].getValue(col);
		}

		avg = sum / rhoO;

		return avg;
	}

	/**
	 * according to type of mutation the algorithm use proper formula and check
	 * boundary for both objective part and strategy part
	 */
	public void mutation() {
		double tu = 1.0 / Math.sqrt(2 * Math.sqrt(dimen));
		double tug = 1.0 / Math.sqrt(2 * dimen);
		double B = 0.0873;
		if (muType == MutationType.UNCORRELATED_TYPE1) {
			for (int i = 0; i < lamdaSize; i++) {
				if (Utils.getRandomValue() < Pm) {
					// strategy
					double sigm = lamdaPop.getIndividual(i).getSigma() * Math.exp(tu * Utils.unformrealRandom());
					lamdaPop.getIndividual(i).setSigma(sigm);
					if (Math.abs(lamdaPop.getIndividual(i).getSigma()) < epsilon)
						lamdaPop.getIndividual(i).setSigma(epsilon);

					// Object Vars
					for (int k = 0; k < dimen; k++) {
						lamdaPop.getIndividual(i).setValue(k, lamdaPop.getIndividual(i).getValue(k)
								+ lamdaPop.getIndividual(i).getSigma() * Utils.unformrealRandom());

						if (lamdaPop.getIndividual(i).getValue(k) < minVal)
							lamdaPop.getIndividual(i).setValue(k, minVal);

						if (lamdaPop.getIndividual(i).getValue(k) > maxVal)
							lamdaPop.getIndividual(i).setValue(k, maxVal);
					}
				}
			}

		} else if (muType == MutationType.UNCORRELATED_TYPE2) {
			for (int i = 0; i < lamdaSize; i++) {
				if (Utils.getRandomValue() < Pm) {
					double randG = Utils.unformrealRandom();
					for (int j = 0; j < dimen; j++) {
						// strategy part
						lamdaPop.getIndividual(i).setSigmaAt(j, lamdaPop.getIndividual(i).getSigmaAt(j)
								* Math.exp(tug * randG + tu * Utils.unformrealRandom()));

						if (Math.abs(lamdaPop.getIndividual(i).getSigmaAt(j)) < epsilon)
							lamdaPop.getIndividual(i).setSigmaAt(j, epsilon);

						// objective part
						lamdaPop.getIndividual(i).setValue(j, lamdaPop.getIndividual(i).getValue(j)
								+ lamdaPop.getIndividual(i).getSigmaAt(j) * Utils.unformrealRandom());

						if (lamdaPop.getIndividual(i).getValue(j) < minVal)
							lamdaPop.getIndividual(i).setValue(j, minVal);

						if (lamdaPop.getIndividual(i).getValue(j) > maxVal)
							lamdaPop.getIndividual(i).setValue(j, maxVal);
					}
				}
			}
		} else if (muType == MutationType.CORRELATED) {

			// Sigmas Mutatios
			for (int i = 0; i < lamdaSize; i++) {
				if (Utils.getRandomValue() < Pm) {
					double randG = Utils.unformrealRandom();
					for (int j = 0; j < dimen; j++) {
						lamdaPop.getIndividual(i).setSigmaAt(j, lamdaPop.getIndividual(i).getSigmaAt(j)
								* Math.exp(tug * randG + tu * Utils.unformrealRandom()));

						if (Math.abs(lamdaPop.getIndividual(i).getSigmaAt(j)) < epsilon)
							lamdaPop.getIndividual(i).setSigmaAt(j, epsilon);
					}
				}
			}
			double tmpLen = dimen * (dimen - 1) / 2;
			// alphas Mutatios
			for (int i = 0; i < lamdaSize; i++) {
				if (Utils.getRandomValue() < Pm) {
					for (int j = 0; j < tmpLen; j++) {
						lamdaPop.getIndividual(i).setAngleAt(j,
								lamdaPop.getIndividual(j).getAngleAt(j) + B * Utils.unformrealRandom());
					}
				}
			}
			// Vars Mutations
			for (int i = 0; i < lamdaSize; i++) {
				if (Utils.getRandomValue() < Pm) {
					Individual tmpInd = lamdaPop.getIndividual(i);
					double[] noise = getNoiceNDime(tmpInd);
					for (int j = 0; j < dimen; j++) {
						double value = tmpInd.getValue(j) + noise[j];
						if (value < minVal)
							value = minVal;

						if (value > maxVal)
							value = maxVal;

						tmpInd.setValue(j, value);
					}
				}
			}
		}
	}

	/**
	 * get a n-dimen noice with rotate approach
	 * 
	 * @param target target individual
	 * @return noices
	 */
	public double[] getNoiceNDime(Individual target) {
		double[] N = new double[dimen];

		for (int i = 0; i < dimen; i++) {
			N[i] = target.getSigmaAt(i) * Utils.unformrealRandom();
		}

		int Nq = dimen * (dimen - 1) / 2;
		Nq--;
		int n1 = 0, n2 = 0, n = dimen - 1;
		for (int k = 0; k < n; k++) {
			n1 = n - k;
			n2 = n;
			for (int i = 0; i < k; i++) {
				double tmpAngle = target.getAngleAt(Nq);
				N[n1] = N[n1] * Math.cos(tmpAngle) - N[n2] * Math.sin(tmpAngle);
				N[n2] = N[n1] * Math.sin(tmpAngle) - N[n2] * Math.cos(tmpAngle);
				n2--;
				Nq--;
			}
		}

		return N;
	}

	/**
	 * there are two type of selection in this algorithm Comma and Plus
	 * 
	 * @param sType type of selection
	 */
	public void Selection(SelectionType sType,int genNumber) {

		if (sType == SelectionType.COMMA) {

			for (int i = 0; i < lamdaPop.size(); i++) {
				lamdaPop.getIndividual(i).setFitness(getFitness(lamdaPop.getIndividual(i).getVars()));
			}

			lamdaPop.sortPop();

			for (int i = 0; i < muSize; i++) {
				muPop.setIndividual(i, lamdaPop.getIndividual(i));
			}

		} else {
			Population lamdMuPop = new Population(muSize + lamdaSize);
			for (int i = 0; i < muSize; i++) {
				lamdMuPop.setIndividual(i, muPop.getIndividual(i));
			}
			for (int j = 0; j < lamdaSize; j++) {
				lamdMuPop.setIndividual(j + muSize, lamdaPop.getIndividual(j));
			}

			for (int i = 0; i < lamdMuPop.size(); i++) {
				lamdMuPop.getIndividual(i).setFitness(getFitness(lamdMuPop.getIndividual(i).getVars()));
				;
			}

			lamdMuPop.sortPop();

			for (int i = 0; i < muSize; i++) {
				muPop.setIndividual(i, lamdMuPop.getIndividual(i));
			}
		}
		if (muPop.getIndividual(0).getFitness() < gbest.getFitness()) {
			gbest.setFitness(muPop.getIndividual(0).getFitness());
			gbest.setVars(muPop.getIndividual(0).getVars());
			gbest.setGeneration(genNumber);
		}
	}

	/**
	 * eval function computation for an individual
	 * 
	 * @param vars objective part of individuals
	 * @return fitness value
	 */
	public double getFitness(double[] vars) {
		double fitness = -1;
		switch (functionEval) {
		case ACKLEY:
			fitness = OptimizationTestFunction.Ackley(vars);
			break;
		case RASTRIGIN:
			fitness = OptimizationTestFunction.Rastrigin(vars);
			break;
		case ROSENBROCK:
			fitness = OptimizationTestFunction.Rosenbrock(vars);
			break;
		case SCHWEFE:
			fitness = OptimizationTestFunction.Schwefel(vars);
			break;
		case SPHERE:
			fitness = OptimizationTestFunction.Sphere(vars);
			break;
		}

		return fitness;
	}

	public static class GBest {
		private double fitness = Double.MAX_VALUE;
		private double[] vars;
		private int generation;

		public double getFitness() {
			return fitness;
		}

		public void setFitness(double fitness) {
			this.fitness = fitness;
		}

		public double[] getVars() {
			return vars;
		}

		public void setVars(double[] vars) {
			this.vars = vars;
		}

		public int getGeneration() {
			return generation;
		}

		public void setGeneration(int generation) {
			this.generation = generation;
		}

		@Override
		public String toString() {
			return String.format("Best Fitness in generation %d is : %.8f%nTraget individual(Objective Part) %s%n",generation, fitness,
					Arrays.toString(vars));
		}
	}

}
