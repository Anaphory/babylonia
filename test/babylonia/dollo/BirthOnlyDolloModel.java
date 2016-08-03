package babylonia.dollo;

import java.util.Arrays;
import java.util.Collection;

import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

@RunWith(Parameterized.class)
public class BirthOnlyDolloModel extends DolloModelTest {
	@Parameters
	public static Collection<Object[]> data() {
		return Arrays.asList(new Object[][] { { 0, 0, 0.0, 0.0 }, { 0, 1, 0.0, 0.0 }, { 1, 0, 0.0, 0.0 },
				{ 1, 1, 0.0, 1.0 }, { 0, 0, 1.0, 0.0 }, { 0, 1, 1.0, 0.0 }, { 1, 0, 1.0, 0.0 }, { 1, 1, 1.0, 1.0 } });
	}

	public BirthOnlyDolloModel(Integer observation1, Integer observation2, Double aliveInEquilibrium,
			Double likelihood) {
		super(observation1, observation2, aliveInEquilibrium, 100.0, 1e-11, likelihood);
	}

}
