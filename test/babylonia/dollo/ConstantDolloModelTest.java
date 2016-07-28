package babylonia.dollo;

import java.util.Arrays;
import java.util.Collection;

public class ConstantDolloModelTest extends DolloModelTest {
	public static Collection<Object[]> data() {
		return Arrays.asList(new Object[][] { { 0, 0, 0.0, 1.0 }, { 0, 1, 0.0, 0.0 }, { 1, 0, 0.0, 0.0 },
				{ 1, 1, 0.0, 0.0 }, { 0, 0, 1.0, 0.0 }, { 0, 1, 1.0, 0.0 }, { 1, 0, 1.0, 0.0 }, { 1, 1, 1.0, 1.0 } });
	}

	public ConstantDolloModelTest(Integer observation1, Integer observation2, Double aliveInEquilibrium,
			Double likelihood) {
		super(observation1, observation2, aliveInEquilibrium, likelihood);
	}

}
