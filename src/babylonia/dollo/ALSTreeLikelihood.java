/*
 * ALSTreeLikelihood.java
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


import beast.core.Description;
import beast.core.Input;


@Description("Treelikelihood for running the Multi-State Stochastic Dollo process")
public class ALSTreeLikelihood extends AbstractObservationProcess {
    public Input<AbstractObservationProcess> opInput = new Input<AbstractObservationProcess>("observationprocess", "description here");

    protected AbstractObservationProcess observationProcess;

    @Override
    public void initAndValidate() {
        observationProcess = opInput.get();
        // ensure TreeLikelihood initialises the partials for tips
        m_useAmbiguities.setValue(true, this);
        super.initAndValidate();
    }

	@Override
	public double calculateLogTreeWeight() {
		// TODO Auto-generated method stub
		return opInput.get().calculateLogTreeWeight();
	}

	@Override
	void setNodePatternInclusion() {
		// TODO Auto-generated method stub
		opInput.get().setNodePatternInclusion();
	}
}
