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
public class ConstantDolloModelTest {
	@Parameters
	public static Collection<Object[]> data() {
		return Arrays.asList(
				new Object[][] { { 0, 0, 0.0, 1.0 }, { 0, 1, 0.0, 0.0 }, { 1, 0, 0.0, 0.0 }, { 1, 1, 0.0, 0.0 },
						{ 0, 0, 1.0, 0.0 }, { 0, 1, 1.0, 0.0 }, { 1, 0, 1.0, 0.0 }, { 1, 1, 1.0, 1.0 } });
	}

	private int[] observations;
	private Double expectedLikelihood;
	protected AnyTipObservationProcessLikelihood dollo;

	public ConstantDolloModelTest(Integer observation1, Integer observation2, Double aliveInEquilibrium,
			Double likelihood) {
		Tree tree = new TreeParser("(A:1,B:1):1");
		observations = new int[] { observation1, observation2 };
		Sequence s1 = new Sequence("A", String.valueOf(observations[0]));
		Sequence s2 = new Sequence("B", String.valueOf(observations[1]));
		Alignment alignment = new Alignment(Arrays.asList(new Sequence[] { s1, s2 }), "binary");

		RealParameter zero = new RealParameter(new Double[] { 1e-11 });
		RealParameter one = new RealParameter(new Double[] { 1.0 });

		dollo = new AnyTipObservationProcessLikelihood();
		SiteModel sites = new SiteModel();
		MutationDeathModel subst = new MutationDeathModel();
		Frequencies freq = new Frequencies(new Double[] { 1.0 - aliveInEquilibrium, aliveInEquilibrium });
		subst.initByName("frequencies", freq, "deathprob", zero);
		sites.initByName("shape", "1.0", "substModel", subst);
		dollo.initByName("tree", tree, "data", alignment, "siteModel", sites, "branchRateModel", new StrictClockModel(),
				"mu", zero, "lambda", one, "integrateGainRate", true);

		expectedLikelihood = likelihood;
	}

	@Test
	public void testTrivialNodePartials() {

		for (int v = 0; v < 2; ++v) {
			for (int i = 0; i < observations.length; ++i) {
				assertEquals((observations[i] == v) ? 1.0 : 0.0, dollo.getNodePartials(i)[v], 1e-7);
			}
		}
	}

	@Test
	public void testCalculateLogP() {
		assertEquals(expectedLikelihood, Math.exp(dollo.calculateLogP()), 1e-8);
	}
}
