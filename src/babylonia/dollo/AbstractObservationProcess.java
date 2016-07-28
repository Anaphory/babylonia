/*
 * AbstractObservationProcess.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond,
 * Andrew Rambaut, Marc Suchard and Alexander V. Alekseyenko
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package babylonia.dollo;

import java.util.HashSet;
import java.util.Set;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.math.GammaFunction;

@Description("Abstract Observation Process defines how the integration of gain events is done along the tree."
		+ "Specific instances should define how the data is collected.")
@Citation("Alekseyenko, AV., Lee, CJ., Suchard, MA. Wagner and Dollo: a stochastic duet"
		+ "by composing two parsimonious solos. Systematic Biology 2008 57(5): 772 - 784; doi:"
		+ "10.1080/10635150802434394. PMID: 18853363")
abstract public class AbstractObservationProcess extends TreeLikelihood {
	public Input<RealParameter> muInput = new Input<RealParameter>("mu",
			"instantaneous per capita loss rate of the character", Validate.REQUIRED);
	public Input<RealParameter> lamInput = new Input<RealParameter>("lam", "instantaneous rate of the Poisson process"
			+ " over the evolutionary tree realizing the gain of characters");
	public Input<Boolean> integrateGainRateInputInput = new Input<Boolean>("integrateGainRate", "description here",
			false);

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub

	}

	protected boolean[] nodePatternInclusion;
	protected boolean[] storedNodePatternInclusion;
	protected double[] cumLike;
	protected double[] nodePartials;
	protected double[] nodeLikelihoods;
	protected int nodeCount;
	protected int patternCount;
	protected int stateCount;
	protected TreeInterface treeModel;
	protected Alignment patterns;
	protected int[] patternWeights;
	protected RealParameter mu;
	protected RealParameter lam;

	// update control variables
	protected boolean weightKnown;
	protected double logTreeWeight;
	protected double storedLogTreeWeight;
	private double gammaNorm;
	private double totalPatterns;
	protected int deathState;
	protected SiteModel siteModel;
	private double logN;
	protected boolean nodePatternInclusionKnown = false;
	BranchRateModel branchRateModel;

	public void initAndValidate(String Name, TreeInterface treeModel, Alignment patterns, SiteModelInterface siteModel,
			BranchRateModel branchRateModel, RealParameter mu, RealParameter lam, boolean integrateGainRate) {
		// super(Name);
		this.treeModel = treeModel;
		this.patterns = patterns;
		this.mu = mu;
		this.lam = lam;
		this.siteModel = (SiteModel) siteModel;
		if (branchRateModel != null) {
			this.branchRateModel = branchRateModel;
		} else {
			this.branchRateModel = new StrictClockModel();
		}
		// addModel(treeModel);
		// addModel(siteModel);
		// addModel(this.branchRateModel);
		// addVariable(mu);
		// addVariable(lam);

		nodeCount = treeModel.getNodeCount();
		stateCount = patterns.getDataType().getStateCount();
		this.patterns = patterns;
		patternCount = patterns.getPatternCount();
		patternWeights = patterns.getWeights();
		totalPatterns = 0;
		for (int i = 0; i < patternCount; ++i) {
			totalPatterns += patternWeights[i];
		}
		logN = Math.log(totalPatterns);

		gammaNorm = -GammaFunction.lnGamma(totalPatterns + 1);

		try {
			this.deathState = ((MutationDeathType) patterns.getDataType()).DEATHSTATE;
		} catch (ClassCastException e) {
			this.deathState = 0;
		}

		setNodePatternInclusion();
		cumLike = new double[patternCount];
		nodeLikelihoods = new double[patternCount];
		weightKnown = false;

		this.integrateGainRate = integrateGainRate;
	}

	public RealParameter getMuParameter() {
		return mu;
	}

	public RealParameter getLamParameter() {
		return lam;
	}

	private double calculateSiteLogLikelihood(int site, double[] partials, double[] frequencies) {
		int v = site * stateCount;
		double sum = 0.0;
		for (int i = 0; i < stateCount; i++) {
			sum += frequencies[i] * partials[v + i];
		}
		return Math.log(sum);
	}

	public final double nodePatternLikelihood(double[] freqs, ALSTreeLikelihood likelihoodCore) {
		int i, j;
		double logL = this.gammaNorm;

		double birthRate = this.lam.getValue(0);
		double logProb;
		if (!this.nodePatternInclusionKnown)
			this.setNodePatternInclusion();
		if (this.nodePartials == null) {
			this.nodePartials = new double[this.patternCount * this.stateCount];
		}

		double averageRate = this.getAverageRate();

		for (j = 0; j < patternCount; ++j)
			this.cumLike[j] = 0;

		for (i = 0; i < nodeCount; ++i) {
			// get partials for node i
			likelihoodCore.getNodePartials(i, this.nodePartials);
			/*
			 * multiply the partials by equilibrium probs this part could be
			 * optimized by first summing and then multiplying by equilibrium
			 * probs
			 */
			// likelihoodCore.calculateLogLikelihoods(nodePartials, freqs,
			// nodeLikelihoods); // MAS Removed
			logProb = Math.log(this.getNodeSurvivalProbability(i, averageRate));

			for (j = 0; j < this.patternCount; ++j) {
				if (this.nodePatternInclusion[i * patternCount + j]) {
					// cumLike[j] += Math.exp(nodeLikelihoods[j] + logProb); //
					// MAS Replaced with line below
					cumLike[j] += Math.exp(this.calculateSiteLogLikelihood(j, this.nodePartials, freqs) + logProb);
				}
			}
		}

		double ascertainmentCorrection = this.getAscertainmentCorrection(this.cumLike);
		// System.err.println("AscertainmentCorrection:
		// "+ascertainmentCorrection);

		for (j = 0; j < patternCount; ++j) {
			logL += Math.log(this.cumLike[j] / ascertainmentCorrection) * this.patternWeights[j];
		}

		double deathRate = this.mu.getValue(0);

		double logTreeWeight = this.getLogTreeWeight();

		if (integrateGainRate) {
			logL -= gammaNorm + logN + Math.log(-logTreeWeight * deathRate / birthRate) * this.totalPatterns;
		} else {
			logL += logTreeWeight + Math.log(birthRate / deathRate) * this.totalPatterns;
		}
		return logL;
	}

	protected double getAscertainmentCorrection(double[] patternProbs) {
		double excludeProb = 0;
		Set<Integer> excludeIndices = new HashSet<Integer>();
	    if (patterns.isAscertained) {
			excludeIndices = patterns.getExcludedPatternIndices();
		}
		for (int index : excludeIndices) {
			excludeProb += patternProbs[index];
		}
		return 1 - excludeProb;
	}

	final public double getLogTreeWeight() {
		if (!this.weightKnown) {
			this.logTreeWeight = this.calculateLogTreeWeight();
			this.weightKnown = true;
		}

		return this.logTreeWeight;
	}

	abstract public double calculateLogTreeWeight();

	abstract void setNodePatternInclusion();

	final public double getAverageRate() {
		if (!averageRateKnown) {
			double avgRate = 0.0;
			double proportions[] = siteModel.getCategoryProportions(null);
			for (int i = 0; i < siteModel.getCategoryCount(); ++i) {
				avgRate += proportions[i] * siteModel.getRateForCategory(i, null);
			}
			averageRate = avgRate;
			averageRateKnown = true;
		}
		return averageRate;
	}

	public double getNodeSurvivalProbability(int index, double averageRate) {
		Node node = treeModel.getNode(index);
		Node parent = node.getParent();

		if (parent == null)
			return 1.0;

		final double deathRate = mu.getValue(0) * averageRate; // getAverageRate();
		final double branchRate = branchRateModel.getRateForBranch(node);
		// Get the operational time of the branch
		final double branchTime = branchRate * node.getLength();
		return 1.0 - Math.exp(-deathRate * branchTime);
	}

	@Override
	public boolean requiresRecalculation() {
		if (mu.somethingIsDirty()) {
			averageRateKnown = false;
			weightKnown = false;
			nodePatternInclusionKnown = false;
		}
		if (siteModel.isDirtyCalculation()) {
			averageRateKnown = false;
		}
		if (treeModel.somethingIsDirty()) {
			weightKnown = false;
			nodePatternInclusionKnown = false;
		}
		return true;
	}

	// protected void handleModelChangedEvent(Model model, Object object, int
	// index) {
	// if (model == siteModel) {
	// averageRateKnown = false;
	// }
	// if (model == treeModel || model == siteModel || model == branchRateModel)
	// {
	// weightKnown = false;
	// }
	// if (model == treeModel) {
	// if (object instanceof TreeModel.TreeChangedEvent) {
	// if (((TreeModel.TreeChangedEvent) object).isTreeChanged()) {
	// nodePatternInclusionKnown = false;
	// }
	// }
	// }
	// }
	//
	// protected final void handleVariableChangedEvent(Variable variable, int
	// index, Parameter.ChangeType type) {
	// if (variable == mu || variable == lam) {
	// weightKnown = false;
	// } else {
	// System.err.println("AbstractObservationProcess: Got unexpected parameter
	// changed event. (Parameter = " + variable + ")");
	// }
	// }
	//
	@Override
	public void store() {
		// storedAverageRate = averageRate;
		storedLogTreeWeight = logTreeWeight;
		System.arraycopy(nodePatternInclusion, 0, storedNodePatternInclusion, 0, storedNodePatternInclusion.length);
	}

	@Override
	public void restore() {
		// averageRate = storedAverageRate;
		averageRateKnown = false;
		logTreeWeight = storedLogTreeWeight;
		boolean[] tmp = storedNodePatternInclusion;
		storedNodePatternInclusion = nodePatternInclusion;
		nodePatternInclusion = tmp;
	}

	protected void acceptState() {
	}

	public void setIntegrateGainRate(boolean integrateGainRate) {
		this.integrateGainRate = integrateGainRate;
	}

	private boolean integrateGainRate = false;

	private double averageRate;
	private boolean averageRateKnown = false;

}