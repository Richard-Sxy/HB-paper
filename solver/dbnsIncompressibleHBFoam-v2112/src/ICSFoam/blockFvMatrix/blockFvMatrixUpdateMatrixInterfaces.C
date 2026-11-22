/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa

    Copyright (C) 2022 Stefano Oliani
-------------------------------------------------------------------------------
License
    This file is part of ICSFOAM.

    ICSFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ICSFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ICSFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "blockFvMatrix.H"
#include "LduInterfaceField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class sourceType, class blockType>
template<class psiType>
void Foam::blockFvMatrix<sourceType, blockType>::initProcessorInterfaces
(
	const Field<psiType>& psiif,
	Field<psiType>& result,
	const LduInterfaceFieldPtrsList<psiType>& interfaces
) const
{
	scalarField interfaceDummyCoeffs(1);

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(interfaces, interfacei)
        {
            if (interfaces.set(interfacei))
            {
                interfaces[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    true,
					this->lduAddr(),
					interfacei,
                    psiif,
					interfaceDummyCoeffs,
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfacei=patchSchedule.size()/2;
            interfacei<interfaces.size();
            interfacei++
        )
        {
            if (interfaces.set(interfacei))
            {
                interfaces[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    true,
					this->lduAddr(),
					interfacei,
                    psiif,
					interfaceDummyCoeffs,
                    Pstream::commsTypes::blocking
                );
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


template<class sourceType, class blockType>
template<class psiType>
void Foam::blockFvMatrix<sourceType, blockType>::updateProcessorInterfaces
(
	Field<sourceType>& result,
	const Pstream::commsTypes commsType,
	const processorFvPatchField<psiType>& procField,
	const fvMesh& mesh,
	const label patchi
) const
{
    if (procField.updatedMatrix()) return;

    // Stock v2112 exposes communication completion via updateInterfaceMatrix.
    // Use one temporary equation per FACE so distinct faces sharing an owner
    // cell are not summed before multiplication by the tensor/vector blocks.
    class faceAddressing final : public lduAddressing
    {
        labelList faces_, empty_;
        lduSchedule schedule_;
    public:
        explicit faceAddressing(const label n) : lduAddressing(n), faces_(n)
        {
            forAll(faces_, i) faces_[i] = i;
        }
        const labelUList& lowerAddr() const override { return empty_; }
        const labelUList& upperAddr() const override { return empty_; }
        const labelUList& patchAddr(const label) const override { return faces_; }
        const lduSchedule& patchSchedule() const override { return schedule_; }
    };

    const faceAddressing addressing(procField.size());
    Field<psiType> neighbour(procField.size(), Zero);
    const scalarField unitCoeffs(procField.size(), 1.0);
    const Field<psiType> unusedInternal;
    // Native implementation scatters with !add: false yields +1*neighbour.
    // It waits for the original requests, transforms once, and sets updatedMatrix.
    procField.updateInterfaceMatrix
    (
        neighbour, false, addressing, 0, unusedInternal, unitCoeffs, commsType
    );

    if (this->interfacesUpper().size())
    {
        const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
        const Field<blockType>& coeffs = this->interfacesUpper()[patchi];
        forAll(neighbour, facei)
            result[faceCells[facei]] += dot(coeffs[facei], neighbour[facei]);
    }

}
// ************************************************************************* //
