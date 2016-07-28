package babylonia.dollo;

import java.util.ArrayList;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.datatype.DataType;

/**
 * Package: MutationDeathType Description: Time: 1:09:40 PM
 */
@Description("Data type for mutation death models including Multi-State Stochastic Dollo")
public class MutationDeathType extends DataType.Base {
	public Input<String> deathCharInput = new Input<String>("deathChar",
			"character representing death state (default 0)", "0");
	public Input<DataType.Base> dataTypeInput = new Input<DataType.Base>("dataType",
			"base datatype, extended by death char");
	public Input<String> extantCodeInput = new Input<String>("extantCode",
			"character representing live state if no existing datatype is extended", Validate.XOR, dataTypeInput);
	// FIXME Why deathChar but extantCode?
	protected static String DESCRIPTION = "MutationDeathType";

	public int DEATHSTATE = 0;

	@Override
	public void initAndValidate() {
		// FIXME This is VERY naive
		char deathCode = deathCharInput.get().charAt(0);
		if (extantCodeInput.get() != null) {
			char extantCode = extantCodeInput.get().charAt(0);

			int[][] x = { { 0 }, // 1
					{ 1 }, // 0
					{ 0, 1 }, // -
					{ 0, 1 }, // ?
			};
			stateCount = 2;
			mapCodeToStateSet = x;
			codeLength = 1;
			codeMap = "" + extantCode + deathCode + GAP_CHAR + MISSING_CHAR;
			DEATHSTATE = 1;
		} else {
			DataType.Base dataType = dataTypeInput.get();
			if (dataType.getCodeLength() != 1) {
				throw new IllegalArgumentException("MutationDeathType only works with code length 1 data types");
			}
			codeLength = 1;
			if (dataType.getCodeMap().contains(String.valueOf(deathCode))) {
				throw new IllegalArgumentException(
						"Death code " + deathCode + " is already a valid code in data type " + dataType.getCodeMap());
			}

			stateCount = dataType.getStateCount() + 1;
			ArrayList<int[]> mapCodeToStateSetList = new ArrayList<int[]>(dataType.getStateCount());
			DEATHSTATE = 0;
			mapCodeToStateSetList.add(new int[] { stateCount - 1 });
			int i = 0;
			try {
				while (true) {
					mapCodeToStateSetList.add(dataType.getStatesForCode(i));
					++i;
				}
			} catch (IndexOutOfBoundsException e) {
				// Do nothing: This exception means the end of the copy process
			}
			mapCodeToStateSetList.toArray(mapCodeToStateSet);
			codeMap = "" + deathCode + dataType.getCodeMap();
			DEATHSTATE = stateCount - 1;
		}
	}

	@Override
	public String getTypeDescription() {
		return "MutationDeathType";
		// FIXME Should this not use the DESCRIPTION?
		// FIXME Should this not refer to the underlying DataType?
	}

}
