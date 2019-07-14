package es;

public class Individual {
	private int dimen;
	private double[] vars;
	private double sigma;
	private double[] sigmaVec = null;
	private double fitness = Double.NaN;
	private double[] alpha = null;

	public Individual(int dimen, double minValue, double maxValue, MutationType type) {
		this.dimen = dimen;
		vars = new double[dimen];
		for (int i = 0; i < dimen; i++) {
			vars[i] = Utils.unformIntRandom(maxValue, minValue);
		}
		if (type == MutationType.UNCORRELATED_TYPE1) {
			sigma = Utils.unformrealRandom();
		} else if (type == MutationType.UNCORRELATED_TYPE2) {
			sigmaVec = new double[dimen];
			for (int i = 0; i < dimen; i++) {
				sigmaVec[i] = Utils.unformrealRandom();
			}
		} else if (type == MutationType.CORRELATED) {
			sigmaVec = new double[dimen];
			for (int i = 0; i < dimen; i++) {
				sigmaVec[i] = Utils.unformrealRandom();
			}
			alpha = new double[(dimen * (dimen - 1)) / 2];
			for (int i = 0; i < alpha.length; i++) {
				alpha[i] = Utils.unformIntRandom(Math.PI, -Math.PI);
			}
		}

	}

	public void setAngleAt(int index, double angle) {
		alpha[index] = getValidAngle(angle);
	}

	public double getAngleAt(int index) {
		return alpha[index];
	}

	public double getValidAngle(double angle) {
		double value = Math.abs(angle);
		if (Double.compare(value, Math.PI) > 0) {
			value = angle - 2 * Math.PI * Math.signum(angle);
			return value;
		}

		return angle;

	}

	public void setSigmaVec(double[] sigmaVec) {
		this.sigmaVec = sigmaVec;
	}

	public double getSigmaAt(int index) {
		return sigmaVec[index];
	}

	public void setSigmaAt(int index, double value) {
		sigmaVec[index] = value;
	}

	public double[] getSigmaVec() {
		return sigmaVec;
	}

	public void setFitness(double fitness) {
		this.fitness = fitness;
	}

	public double getFitness() {
		return fitness;
	}

	public void setValue(int index, double value) {
		this.vars[index] = value;
	}

	public double getValue(int index) {
		return vars[index];
	}

	public int getDimen() {
		return dimen;
	}

	public double[] getVars() {
		return vars;
	}

	public void setVars(double[] vars) {
		this.vars = vars;
	}

	public double getSigma() {
		return sigma;
	}

	public void setSigma(double sigma) {
		this.sigma = sigma;
	}

	@Override
	public String toString() {
		StringBuilder str = new StringBuilder();
		str.append("< ");
		for (double item : vars) {
			str.append(item + ",");
		}
		if (alpha != null) {
			str.append("sigma:" + sigmaVec[0]);
			for (int i = 1; i < sigmaVec.length; i++)
				str.append("," + sigmaVec[i]);
			str.append(",Alpha: " + alpha[0]);
			for (int i = 1; i < alpha.length; i++)
				str.append("," + alpha[i]);

			str.append(" >");
		} else if (sigmaVec != null) {
			str.append("sigma:" + sigmaVec[0]);
			for (int i = 1; i < sigmaVec.length; i++)
				str.append("," + sigmaVec[i]);

			str.append(" >");
		} else {
			str.append("sigma:" + sigma + " > ");
		}

		return str.toString();
	}
}
