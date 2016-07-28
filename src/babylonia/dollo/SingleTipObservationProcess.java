/*
 * SingleTipObservationProcess.java
 *
 * Copyright (C) 2002-2016 Gereon Kaiping, Alexei Drummond,
 * Andrew Rambaut, Marc Suchard and Alexander V. Alekseyenko
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

import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Taxon;

@Description("Observation process for Multi-State Stochastic Dollo model. Defines a data collection process where the traits must be present in a specific tip node.")
public class SingleTipObservationProcess extends AnyTipObservationProcess {

	public Input<Taxon> theTip = new Input<Taxon>("taxon", "A taxon in which the traits must be present",
			Input.Validate.REQUIRED);
	Taxon sourceTaxon;

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		this.sourceTaxon = theTip.get();
		// FIXME: Assert that the characters are indeed present in the tip.
	}

	@Override
	public double calculateLogTreeWeight() {
		return -lam.getValue(0) / (getAverageRate() * mu.getValue(0));
	}
}
