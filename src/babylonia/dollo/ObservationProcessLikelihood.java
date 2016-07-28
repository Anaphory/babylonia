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
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.math.GammaFunction;

@Description("Tree Likelihood using an "
		+ "Abstract Observation Process, which defines how the integration of gain events is done along the tree."
		+ "Specific instances should define how the data is collected.")
@Citation("Alekseyenko, AV., Lee, CJ., Suchard, MA. Wagner and Dollo: a stochastic duet"
		+ "by composing two parsimonious solos. Systematic Biology 2008 57(5): 772 - 784; doi:"
		+ "10.1080/10635150802434394. PMID: 18853363")
public abstract class ObservationProcessLikelihood extends TreeLikelihood {
	public abstract double calculateLogTreeWeight();
	public abstract void setNodePatternInclusion();

	public Input<RealParameter> muInput = new Input<RealParameter>("mu",
			"instantaneous per capita loss rate of the character", Validate.REQUIRED);
	public Input<RealParameter> lamInput = new Input<RealParameter>("lambda",
			"instantaneous rate of the Poisson process"
					+ " over the evolutionary tree realizing the gain of characters",
			new RealParameter(new Double[] { 1.0 }));


	protected boolean[] nodePatternInclusion;
	private boolean[] storedNodePatternInclusion;
	protected int patternCount;
	protected int stateCount;
	protected TreeInterface tree;
	protected int[] patternWeights;
	protected RealParameter mu;
	protected RealParameter lam;

	// update control variables
	protected boolean weightKnown;
	protected double logTreeWeight;
	protected double storedLogTreeWeight;
	protected double gammaNorm;
	protected double totalPatterns;
	public int deathState;
	protected SiteModel siteModel;
	protected boolean nodePatternInclusionKnown = false;
	protected BranchRateModel branchRateModel;

	@Override
	public void initAndValidate() {
		// ensure TreeLikelihood initialises the partials for tips
		m_useAmbiguities.setValue(true, this);
		super.initAndValidate();
		// FIXME: Explicit cast to SiteModel – Should the methods we call be in
		// the Interface, should GenericTreeLikelihood require a SiteModel, is
		// this a special case for ObservationProcessLikelihood or can we do
		// without these methods?
		initAndValidate("AnyTip", treeInput.get(), dataInput.get(), (SiteModel) siteModelInput.get(),
				branchRateModelInput.get(),
				muInput.get(),
				lamInput.get());
	}

	public void initAndValidate(String Name, TreeInterface treeInterface, Alignment patterns,
			SiteModel siteModelInterface,
			BranchRateModel branchRateModel, RealParameter mu, RealParameter lam) {
		this.tree = treeInterface;
		this.mu = mu;
		this.lam = lam;
		this.siteModel = siteModelInterface;
		if (branchRateModel != null) {
			this.branchRateModel = branchRateModel;
		} else {
			this.branchRateModel = new StrictClockModel();
		}

		stateCount = patterns.getDataType().getStateCount();
		patternCount = patterns.getPatternCount();
		patternWeights = patterns.getWeights();
		totalPatterns = 0;
		for (int i = 0; i < patternCount; ++i) {
			totalPatterns += patternWeights[i];
		}

		gammaNorm = -GammaFunction.lnGamma(totalPatterns + 1);

		try {
			this.deathState = ((MutationDeathType) patterns.getDataType()).DEATHSTATE;
		} catch (ClassCastException e) {
			this.deathState = 0;
		}
		// FIXME We need to actually check that "0" is in the data type

		setNodePatternInclusion();
		weightKnown = false;
	}
	
	@Override
	public double calculateLogP() {
		// Calculate the partial likelihoods
		super.calculateLogP();
		// get the frequency model
		double[] freqs = ((SiteModel.Base) siteModelInput.get()).substModelInput.get().getFrequencies();
		// let the observationProcess handle the rest
		logP = nodePatternLikelihood(freqs);
		return logP;
	}

	public void getNodePartials(int iNode, double[] fPartials) {
		if (requiresRecalculation()) {
			calculateLogTreeWeight();
			calculateLogP();
		}
		if (beagle != null) {
			throw new RuntimeException(this.getClass().getName() + " does not support Beagle yet.");
			// beagle.beagle.getPartials(beagle.partialBufferHelper.getOffsetIndex(iNode),
			// Beagle.NONE, fPartials);
		} else {
			likelihoodCore.getNodePartials(iNode, fPartials);
		}
	}

	public double[] getNodePartials(int iNode) {
		double[] fPartials = new double[patternCount * stateCount];
		getNodePartials(iNode, fPartials);
		return fPartials;
	}


	public RealParameter getMuParameter() {
		return mu;
	}

	public RealParameter getLamParameter() {
		return lam;
	}

	private double calculateSiteLikelihood(int site, double[] partials, double[] frequencies) {
		int v = site * stateCount;
		double sum = 0.0;
		for (int i = 0; i < stateCount; i++) {
			sum += frequencies[i] * partials[v + i];
		}
		return sum;
	}

	public final double nodePatternLikelihood(double[] freqs) {
		int i, j;
		double logL = gammaNorm;

		if (!nodePatternInclusionKnown)
			setNodePatternInclusion();
		
		double[] nodePartials = new double[patternCount * stateCount];
		double averageRate = getAverageRate();

		double[] cumLike = new double[patternCount];

		for (j = 0; j < patternCount; ++j) {
			cumLike[j] = 0;
		}

		for (i = 0; i < tree.getNodeCount(); ++i) {
			likelihoodCore.getNodePartials(i, nodePartials);
			double nodeSurvival = getNodeLossProbability(i, averageRate);
			if (nodeSurvival == 0) {
				continue;
			}

			for (j = 0; j < patternCount; ++j) {
				if (nodePatternInclusion[i * patternCount + j]) {
					cumLike[j] += calculateSiteLikelihood(j, nodePartials, freqs) * nodeSurvival;
				}
			}
		}

		double ascertainmentCorrection = getAscertainmentCorrection(cumLike);
		// System.err.println("AscertainmentCorrection:
		// "+ascertainmentCorrection);

		for (j = 0; j < patternCount; ++j) {
			logL += Math.log(cumLike[j] / ascertainmentCorrection) * patternWeights[j];
		}

		logL += gainRate();

		return logL;
	}

	protected double getAscertainmentCorrection(double[] patternProbs) {
		double excludeProb = 0;
		Set<Integer> excludeIndices = new HashSet<Integer>();
		if (dataInput.get().isAscertained) {
			excludeIndices = dataInput.get().getExcludedPatternIndices();
		}
		for (int index : excludeIndices) {
			excludeProb += patternProbs[index];
		}
		return 1 - excludeProb;
	}

	final public double getLogTreeWeight() {
		if (!weightKnown) {
			logTreeWeight = calculateLogTreeWeight();
			weightKnown = true;
		}

		return logTreeWeight;
	}

	// should override this method

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

	/**
	 * Calculate the probability that a state which is present in the parent is
	 * lost along the branch leading to a node.
	 * 
	 * @param index:
	 *            the index (in `tree`) of the node i of interest
	 * @param averageRate:
	 *            the average rate of evolution
	 * @return P(i=0|i.parent≠0)
	 */
	public double getNodeLossProbability(int index, double averageRate) {
		Node node = tree.getNode(index);
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
		if (tree.somethingIsDirty()) {
			weightKnown = false;
			nodePatternInclusionKnown = false;
		}
		return true;
	}

	@Override
	public void store() {
		// storedAverageRate = averageRate;
		storedLogTreeWeight = logTreeWeight;
		storedNodePatternInclusion = new boolean[nodePatternInclusion.length];
		System.arraycopy(nodePatternInclusion, 0, storedNodePatternInclusion, 0, nodePatternInclusion.length);
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

	private double averageRate;
	private boolean averageRateKnown = false;

	protected double gainRate() {
		 return getLogTreeWeight() + Math.log(lam.getValue(0) / mu.getValue(0)) * totalPatterns;
	}

}
