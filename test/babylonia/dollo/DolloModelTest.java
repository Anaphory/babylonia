/*
 * ConstantDolloModelTest.java
 *
 * Copyright (C) 2016 Gereon Kaiping
 *
 * This file is part of the Beast2 package extension babylonia.
 * See the COPYING file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * This package (just like Beast2) is free software; you can
 * redistribute it and/or modify it under the terms of the GNU
 * Lesser General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option)
 * any later version.
 *
 * The software is distributed in the hope that it will be useful,
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

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.Collection;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.MutationDeathModel;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;

@RunWith(Parameterized.class)
public class DolloModelTest {
	@Parameters
	public static Collection<Object[]> data() {
		return Arrays.asList(new Object[][] { { 0, 0, 0.0, 1.0 }, { 0, 1, 0.0, 0.0 }, { 1, 0, 0.0, 0.0 },
				{ 1, 1, 0.0, 0.0 }, { 0, 0, 1.0, 0.0 }, { 0, 1, 1.0, 0.0 }, { 1, 0, 1.0, 0.0 }, { 1, 1, 1.0, 1.0 } });
	}

	private int[] observations;
	private Double expectedLikelihood;
	protected AnyTipObservationProcess dollo;

	public DolloModelTest(Integer observation1, Integer observation2, Double aliveInEquilibrium,
			Double likelihood) {
		Tree tree = new TreeParser("(A:1,B:1):1");
		observations = new int[] { observation1, observation2 };
		Sequence s1 = new Sequence("A", String.valueOf(observations[0]));
		Sequence s2 = new Sequence("B", String.valueOf(observations[1]));
		Alignment alignment = new Alignment();
		MutationDeathType dtype = new MutationDeathType();
		dtype.initByName("extantCode", "1");
		alignment.initByName("sequence", Arrays.asList(new Sequence[] { s1, s2 }), "userDataType", dtype);

		RealParameter zero = new RealParameter(new Double[] { 1e-11 });
		RealParameter one = new RealParameter(new Double[] { 1.0 });

		dollo = new AnyTipObservationProcess();
		SiteModel sites = new SiteModel();
		MutationDeathModel subst = new MutationDeathModel();
		// NOTE: The encoding of the basic MutationDeatType is "1"→0 and
		// "0"→1, and the frequencies are noted in ENCODING order, not in
		// CHARACTER order!
		Frequencies freq = new Frequencies(new Double[] { aliveInEquilibrium, 1.0 - aliveInEquilibrium });
		subst.initByName("frequencies", freq, "deathprob", zero);
		sites.initByName("shape", "1.0", "substModel", subst);
		dollo.initByName("tree", tree, "data", alignment, "siteModel", sites, "branchRateModel", new StrictClockModel(),
				"mu", zero, "lam", zero, "integrateGainRate", true);

		expectedLikelihood = likelihood;
	}

	@Test
	public void testTrivialNodePartials() {

		for (int v = 0; v < 2; ++v) {
			for (int i = 0; i < observations.length; ++i) {
				double[] fPartials = new double[2];
				dollo.getNodePartials(i, fPartials);
				// NOTE: The encoding of the basic MutationDeatType is "1"→0 and
				// "0"→1, and the frequency partials are noted in ENCODING
				// order, not in CHARACTER order!
				assertEquals((observations[i] == v) ? 1.0 : 0.0, fPartials[1 - v], 1e-7);
			}
		}
	}

	@Test
	public void testCalculateLogP() {
		assertEquals(expectedLikelihood, Math.exp(dollo.calculateLogP()), 1e-8);
	}
}
