/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be Ubeameful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 UbeamA


\*---------------------------------------------------------------------------*/

#include "phaseCoupling.H"
#include "mathematicalConstants.H"
#include "clock.H"
#include "fvMesh.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::phaseCoupling::checkvolumeConnection(const label startCell, const label currentID)
{
    // Set volume ID
    vofID_[startCell] = currentID;

    // Create dynamic list to store all neighbours
    DynamicList<label> cells(0);
    // Get neighbours of startCell
    cells.append(mesh_.cellCells()[startCell]);

    // Start outer loop
    while (cells.size() > 0)
    {
        // Create new dynamic list to store all new neighbours
        DynamicList<label> neighbours(0);

        // Loop over all current neighbours
        forAll(cells, i)
        {
            label cellI = cells[i];
            // Check if neighbour cell is gas and has no ID yet
            if (alpha_[cellI] < alphaLimit_ && mag(vofID_[cellI] - currentID) > 0.1)
            {
                // Set volume ID
                vofID_[cellI] = currentID;

                // add neighbours to neighboursNeighbours list
                neighbours.append(mesh_.cellCells()[cellI]);
            }
        }
        cells.clear();
        cells.transfer(neighbours);
    }
}

Foam::scalar Foam::phaseCoupling::gaussiankernelfunction
(
     scalar x,
     scalar R
)
{
     scalar sum = 1/pow(sqrt(2*Foam::constant::mathematical::pi)*R,3)*exp(-pow(x,2)*0.5);
     return sum;
}

Foam::scalar Foam::phaseCoupling::voidfraction
(
     scalar x,
     scalar R,
     scalar vcell,
     scalar vb
)
{
     scalar v = vb * gaussiankernelfunction(x,R)/vcell;
     return v;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::phaseCoupling::phaseCoupling
(
    const volVectorField& U,
    volScalarField& alpha,
    volScalarField& p,
    volScalarField& beta,
    bubbleCloud& Cloud_
)
:
    mesh_(U.mesh()),
    dict_
    (
        IOobject
        (
            "cloudProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    U_(U),
    beta_(beta),
    alpha_(alpha),
    p_(p),
    cloud_(Cloud_),
    VoFIDactive_(readBool(dict_.subDict("phaseCoupling").lookup("VoFIDactive"))),
    alphaLimit_(readScalar(dict_.subDict("phaseCoupling").lookup("alphaLimit"))),
    startTime_(readScalar(dict_.subDict("phaseCoupling").lookup("startTime"))),
    nInterval_(readLabel(dict_.subDict("phaseCoupling").lookup("nInterval"))),
    vlim_(readScalar(dict_.subDict("phaseCoupling").lookup("vlim"))),
    vofID_
    (
        IOobject
        (
            "vofID",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("vofID", dimless, scalar(0.0)),
        zeroGradientFvPatchScalarField::typeName
    ),
    tNumber_(0),
    pCenter_(0),
    pRadius_(0),
    pGrids_(0),
    pTracks_(0),
    closeVoF_(0),
    oldMom_(0),
    newMom_(0),
    cVolume_(0.0),
    cddRt_(0.0),
    overlap_(mesh_.nCells()),
    innerVolume_(0),
    interfaceVolume_(0),
    effectZone_
    (
        IOobject
        (
            "effectZone",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("effecZone", dimless, scalar(0.0)),
        zeroGradientFvPatchScalarField::typeName
    ),
    momentumSource_
    (
        IOobject
        (
           "momentumSource",
           mesh_.time().timeName(),
           mesh_,
           IOobject::NO_READ,
           IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0,0,0),vector::zero)
    ),
    alphaEqnCouple_(dict_.subDict("phaseCoupling").lookupOrDefault<Switch>("alphaEqnCouple",false)),
    momentumCouple_(dict_.subDict("phaseCoupling").lookupOrDefault<Switch>("momentumCouple",false)),
    Rcriterion_(dict_.subDict("phaseCoupling").lookupOrDefault<Switch>("Rcriterion",true)),
    NleCriterion_(dict_.subDict("phaseCoupling").lookupOrDefault<Switch>("NleCriterion",false)),
    Nle_(dict_.subDict("phaseCoupling").lookupOrDefault<label>("Nle",10)),
    Rthr_(readScalar(dict_.subDict("phaseCoupling").lookup("Rthr"))),
    NelCriterion_(dict_.subDict("phaseCoupling").lookupOrDefault<Switch>("NelCriterion",false)),
    Nel_(dict_.subDict("phaseCoupling").lookupOrDefault<label>("Nel",20))

{
    Info<<"Rcriterion = "<<Rcriterion_;
    if(Rcriterion_)
    {
        Info<<" Rthr_:  "<<Rthr_;
    }
    Info<<endl;
    Info<<"NleCriterion =  "<<NleCriterion_<<endl;
    if(NleCriterion_)
    {
        Info<<" Nle_:  "<<Nle_;
    }
    Info<<endl;
    Info<<"VolumeCirterion = "<<VoFIDactive_ <<endl;
    if(VoFIDactive_)
    {
        Info<<"vlim_:  "<<vlim_;
    }
    Info<<endl;
    Info<<"NelCriterion =  "<<NelCriterion_<<endl;
    
    if(NelCriterion_)
    {
       Info<<" Nel_:  "<<Nel_;
    }
    Info<<endl;

}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::phaseCoupling::update()
{
//   updateBubble();
   updateBeta(beta_);
   updatevolume();
   writeLagrangianBubble();
}
Foam::word Foam::phaseCoupling::relativePos(const vector& center, const scalar& radius, const label& celli)
{
    const vectorList& pts = mesh_.points();
    labelList tempPoints(mesh_.cellPoints()[celli]);
    label  value = 0.0;
    forAll(tempPoints,k)
    {
        if(mag(pts[tempPoints[k]] - center) <= radius)
        {
            value ++;
        }
    }
    if( value==0 && celli == mesh_.findCell(center))
    {
        return "encircle";
    }
    if(value == tempPoints.size())  //all points are in the bubble
    {
        return "inner";
    }
    else if(value != 0)  //some points in the bubble
    {
        return "interface";
    }
    return "outside";
}

void Foam::phaseCoupling::updateBubble(bool remove)
{
    labelList nParticles(Pstream::nProcs(),0); //List index: procssore id; List command: bubble number in corresponding processore
    nParticles[Pstream::myProcNo()] = cloud_.nParcels();
    reduce(nParticles,sumOp<labelList>());

    label Nparticles = cloud_.nParcels();
    reduce(Nparticles,sumOp<label>());

    pCenter_.clear();
    pRadius_.clear();
    pGrids_.clear();
    oldMom_.clear();
    newMom_.clear();
    closeVoF_.clear();
    innerVolume_.clear();
    interfaceVolume_.clear();
    
    pCenter_.resize(Nparticles,vector(0,0,0));
    pRadius_.resize(Nparticles,0.0);
    pGrids_.resize(Nparticles,0);
    closeVoF_.resize(Nparticles,0);
    pTracks_.resize(Nparticles,0);
    oldMom_.resize(Nparticles,vector(0,0,0));
    newMom_.resize(Nparticles,vector(0,0,0));
    innerVolume_.resize(Nparticles,0.0);
    interfaceVolume_.resize(Nparticles,0.0);
    cVolume_ = 0; 
    cddRt_ = 0; 
    label n = 0;
    forAll(nParticles,i)
    {
        label j = Pstream::myProcNo();
        if(i < j)
        {
            n += nParticles[i];
        }
    }
    label procsStart = n;

    for(bubble& bubble : cloud_)
    {
        pCenter_[n] = bubble.position();
        pRadius_[n] = bubble.R();
        vector U0 = bubble.U0();
        vector U1 = bubble.U();
        scalar R = bubble.R();
        scalar mass = cloud_.rhop()*pow(R,3.0)*Foam::constant::mathematical::pi*4.0/3.0;
        cVolume_ += pow(R,3.0)*Foam::constant::mathematical::pi*4.0/3.0;
        cddRt_ += bubble.ddRt()*4*Foam::constant::mathematical::pi*pow(R,2.0) + 8*Foam::constant::mathematical::pi*R*pow(bubble.dRt(),2.0);
        oldMom_[n] = mass * U0; 
        newMom_[n] = mass * U1;
        n ++;       
    }
    reduce(pCenter_,sumOp<vectorList>());
    reduce(pRadius_,sumOp<scalarList>());
    reduce(oldMom_,sumOp<vectorList>());
    reduce(newMom_,sumOp<vectorList>());
    reduce(cVolume_,sumOp<scalar>());
    reduce(cddRt_,sumOp<scalar>());
    overlap_.clear();
    overlap_.resize(mesh_.nCells());   
    
    effectZone_ = 0.0;

    forAll(pRadius_,i)
    {
        vector center = pCenter_[i];
        scalar radius = pRadius_[i];
        labelList  effectCell(0);
        label  grids = 0;

        forAll(effectZone_,celli)
        {
            if( mag(mesh_.C()[celli]-center) <= 5.0*radius)
            {
                effectCell.append(celli);
            }
        }
        forAll(effectCell,j)
        {
            word topo = relativePos(center,radius,effectCell[j]);
            if(topo == "inner")
            {
                overlap_[effectCell[j]].append(i);
                innerVolume_[i] += mesh_.V()[effectCell[j]];
                grids++;
            }
            else if(topo == "interface")
            {
                overlap_[effectCell[j]].append(i);
                interfaceVolume_[i] += mesh_.V()[effectCell[j]];
                grids++;
            }
            else if(topo == "encircle")
            {
                overlap_[effectCell[j]].append(i);
                grids++;
            }
        }
        pGrids_[i] = grids++;
    }
    reduce(pGrids_,sumOp<labelList>());
    reduce(innerVolume_,sumOp<scalarList>());
    reduce(interfaceVolume_,sumOp<scalarList>());

    //check lagrangian bubble if close to VoF
    const volScalarField& alpha1 = mesh_.lookupObject<const volScalarField>("alpha.water");
    forAll(overlap_,celli)
    {
        labelList bubbles = overlap_[celli];
        if(bubbles.size() != 0)
        {
            forAll(bubbles,i)
            {
                word topo = relativePos(pCenter_[bubbles[i]],pRadius_[bubbles[i]],celli);
                if((topo == "interface" || topo == "inner" || topo == "encircle") && alpha1[celli] < 0.6)
                {
                    closeVoF_[bubbles[i]] = 1;
                }
            }
        }
    }
    reduce(closeVoF_,maxOp<labelList>());
    
    //delete bubble: 1.close to VoF; 2.bigger than Rthr; 3.Occupy cells more than Nle
    if(alphaEqnCouple_ && remove)
    {
        DynamicList<label> removeBubbles(0);
        forAll(pRadius_,i)
        {
            if(pRadius_[i] > Rthr_ && Rcriterion_)
            {
                removeBubbles.append(i);
                pTracks_[i] = 1;
                tNumber_ ++;
            }
            else if(closeVoF_[i] == 1)
            {
                removeBubbles.append(i);
                pTracks_[i] = 1;
                tNumber_ ++;
            }
            else if(pGrids_[i] > Nle_ && NleCriterion_)
            {
                removeBubbles.append(i);
                pTracks_[i] = 1;
                tNumber_ ++;
            }
        }
        reduce(pTracks_,maxOp<labelList>());
        label localIndex = 0;
        for(bubble& bubble : cloud_)
        {
            if(removeBubbles.found(localIndex+procsStart))
            {
                cloud_.deleteParticle(bubble);
                //scalar R = bubble.R();
               // scalar ddRt = bubble.ddRt();
               // scalar dRt = bubble.dRt();
               // cVolume_ -= pow(R,3.0)*Foam::constant::mathematical::pi*4.0/3.0;
               // cddRt_ -= ddRt*4.0*Foam::constant::mathematical::pi*pow(R,2.0) + 8*Foam::constant::mathematical::pi*R*pow(dRt,2.0);
            }
            localIndex++;
        }
    }
           
    //  reduce(cVolume_,sumOp<scalar>());
    // reduce(cddRt_,sumOp<scalar>());

    if(Pstream::master())
    {
         Info<<"update lagrangian bubble successful"<<endl;
         Info<<"translation bubble number: "<<tNumber_<<endl;
    }
}
void Foam::phaseCoupling::updateBeta(volScalarField& beta)
{
   
    beta = 1.0;
    momentumSource_ = dimensionedVector("zero", dimensionSet(1,-2,-2,0,0,0,0),vector::zero);
    
        scalarList Vcellk(pRadius_.size(),0.0);
        forAll(overlap_,cellI)
        {
            labelList bubbles = overlap_[cellI];
            if(bubbles.size() != 0)
            {
                forAll(bubbles,i)
                {
                    scalar R = pRadius_[bubbles[i]];
                    vector posP = pCenter_[bubbles[i]];
                    scalar Vp=pow(2.0*R,3)*constant::mathematical::pi/6.0;
                    vector poscell = mesh_.cellCentres()[cellI];
                    scalar x = sqrt(pow((poscell.component(0)-posP.component(0)),2)+pow((poscell.component(1)-posP.component(1)),2)+pow((poscell.component(2)-posP.component(2)),2));
                    if( x <= 1.16*R && pTracks_[bubbles[i]] == 0)
                    {
                        Vcellk[bubbles[i]] += voidfraction(x/R,R,mesh_.cellVolumes()[cellI],Vp) * mesh_.cellVolumes()[cellI];
                    }
                }
            }
        }
        reduce(Vcellk,sumOp<scalarList>());
        forAll(overlap_,cellI)
        {
            labelList bubbles = overlap_[cellI];
            if(bubbles.size() != 0)
            {
                forAll(bubbles,i)
                {
                    scalar R = pRadius_[bubbles[i]];
                    vector posP = pCenter_[bubbles[i]];
                    scalar Vp=pow(2.0*R,3)*constant::mathematical::pi/6.0;
                    vector momChange = (newMom_[bubbles[i]] - oldMom_[bubbles[i]])/cloud_.deltaTL();
                    label cellme=mesh_.findCell(posP);
                    vector poscell = mesh_.cellCentres()[cellI];
                    scalar x = sqrt(pow((poscell.component(0)-posP.component(0)),2)+pow((poscell.component(1)-posP.component(1)),2)+pow((poscell.component(2)-posP.component(2)),2));
                    if((x <= 1.16* R || cellI == cellme) && pTracks_[bubbles[i]] == 0 )
                    {
                        if(Vp < mesh_.cellVolumes()[cellme] || Vcellk[bubbles[i]] == 0)
                        {
                            Vcellk[bubbles[i]] = mesh_.cellVolumes()[cellI];
                            beta[cellI] -= Vp/(Vcellk[bubbles[i]]+VSMALL);
                            if(momentumCouple_)
                            {        
                                     momentumSource_[cellI] += momChange/(mesh_.cellVolumes()[cellI]+VSMALL);
                            }
                        }
                        else
                        {
                            beta[cellI] -= voidfraction(x/R,R,mesh_.cellVolumes()[cellI],Vp)*Vp/(Vcellk[bubbles[i]]+VSMALL);
                            if(momentumCouple_)
                            {
                                  momentumSource_[cellI] += voidfraction(x/R,R,mesh_.cellVolumes()[cellI],Vp)*1.0/(Vcellk[bubbles[i]]+VSMALL)*mesh_.cellVolumes()[cellI] * momChange/(mesh_.cellVolumes()[cellI]+VSMALL);
                            }
                        }
                    }
                }
            }
	}
        beta.boundaryFieldRef() = 1.0;
        momentumSource_.boundaryFieldRef() = vector::zero;
        beta = min(beta,scalar(1.0));
        beta = max(beta,scalar(0.0));
        if(Pstream::master())
        {
            Info<<"update beta and momentum successful"<<endl;
        }

}

void Foam::phaseCoupling::updateAlphaSource(volScalarField& Sa)
{
    const volScalarField& alpha1 = mesh_.lookupObject<const volScalarField>("alpha.water");
    volScalarField alpha_Sa
    (
        IOobject
        (
            "alpha_Sa",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        -alpha1/mesh_.time().deltaT()*mesh_.time().deltaT().value()
    );
 
    Sa = 0.0/mesh_.time().deltaT();
    if(alphaEqnCouple_)
    {
        forAll(overlap_,celli)
        {
            labelList bubbles = overlap_[celli];
            if(bubbles.size() != 0)
            {
                forAll(bubbles,i)
                {
                   if((pRadius_[bubbles[i]] > Rthr_ && Rcriterion_ )|| (pGrids_[bubbles[i]] > Nle_ && NleCriterion_) || closeVoF_[bubbles[i]] == 1)
                    {
                        word topo = relativePos(pCenter_[bubbles[i]],pRadius_[bubbles[i]],celli);
                        scalar bubbleVolume = 4.0/3.0*constant::mathematical::pi*pow(pRadius_[bubbles[i]],3.0);
                        if(topo == "inner")
                        {
                            Sa[celli] -= 1.0;
                        }
                        else if(topo == "encircle")
                        {
                            Sa[celli] -= bubbleVolume/mesh_.V()[celli];
                        }
                        else if(topo == "interface")
                        {
                            Sa[celli] -= (bubbleVolume - innerVolume_[bubbles[i]])/interfaceVolume_[bubbles[i]];
                        }
                    }
                }
            }
        }
        Sa.boundaryFieldRef() = 1.0;
        Sa = min(Sa,dimensionedScalar("maxSa",Sa.dimensions(),0.0));
        Sa = max(Sa,alpha_Sa);
        Sa = Sa/mesh_.time().deltaT().value();
    }
}

const volVectorField&  Foam::phaseCoupling::momentumSource()
{
    return momentumSource_;
}
void Foam::phaseCoupling::initial()
{
    if(Pstream::master())
    {
        fileName outputDir = mesh_.time().path()/"postProcessing";
        mkDir(outputDir);
        outputPtrR_.reset(new OFstream(outputDir/"lagrangianBubble_R.dat"));
        outputPtrX_.reset(new OFstream(outputDir/"lagrangianBubble_X.dat"));
        outputPtrV_.reset(new OFstream(outputDir/"lagrangianBubble_V.dat"));
        outputPtrddRt_.reset(new OFstream(outputDir/"lagrangianBubble_ddRt_.dat"));
        if(Pstream::master())
        {
            Info<<"initial file writer pointer successful"<<endl;
        }
    }
}

void Foam::phaseCoupling::writeLagrangianBubble()
{
    if(pRadius_.size() != 0)
    {
        if(Pstream::master())
        {
            outputPtrR_()<<mesh_.time().timeName()<<",";
            outputPtrX_()<<mesh_.time().timeName()<<",";           
            outputPtrV_()<<mesh_.time().timeName()<<",";           
            outputPtrddRt_()<<mesh_.time().timeName()<<",";           
            forAll(pCenter_,i)
            {
                outputPtrX_()<<pCenter_[i]<<",";        
            }
            outputPtrX_()<<endl;
            forAll(pRadius_,i)
            {
                outputPtrR_()<<pRadius_[i]<<",";
            }
            outputPtrR_()<<endl;    
            outputPtrV_()<<cVolume_<<endl;    
            outputPtrddRt_()<<cddRt_<<endl;    
        }
    }
}

// trannsition with pressure higher than limit pressure
void Foam::phaseCoupling::updatevolume()
{
    // Reset volumeID and damping/source fields
    vofID_ *= 0.0;

    if
    (
        mesh_.time().timeIndex() % nInterval_ != 0
     || mesh_.time().time().value() < startTime_
    )
    {
        return;
    }

    // Global time step
    dimensionedScalar deltaT = mesh_.time().deltaT();

    // Local-to-global cell reference
    globalIndex globalNumbering (mesh_.nCells());
    // Create the volumeID field
    forAll(mesh_.cells(), cellI)
    {
        // Convert local processor cellIDs to global cellIDs
        label globalCellID = globalNumbering.toGlobal(cellI);

        // Cell is considered gas and not yet marked
        if (alpha_[cellI] < alphaLimit_ && mag(vofID_[cellI]) < 0.1)
        {
            // Check surrounding cells
            checkvolumeConnection(cellI, globalCellID+1);
        }
    }
    vofID_.correctBoundaryConditions();


    // Correct for parallel processing
    // Number of corrections
    label corr = Foam::log(scalar(Pstream::nProcs())) / Foam::log(2.0) + 1;
    for(int i=0; i<=corr; i++)
    {
        forAll(mesh_.boundary(), patchI)
        {
            // Current patch
            const fvPatch& pPatch = mesh_.boundary()[patchI];

            // check if patch is coupled, e.g. processor patch
            if (pPatch.coupled())
            {
                // loop over all faces of processor patch
                forAll(pPatch, faceI)
                {
                    // get cell values
                    const scalar ID_own = vofID_.boundaryField()[patchI].patchInternalField()()[faceI];
                    const scalar ID_nei = vofID_.boundaryField()[patchI].patchNeighbourField()()[faceI];

                    // check, if volumeIDs on both sides are vaper
                    if ( (ID_nei > 0.1) && (ID_own > 0.1) )
                    {
                        // check, if volumeIDs on both sides are different
                        if ( mag(ID_nei - ID_own) > 0.1 )
                        {
                            scalar minID = Foam::min(ID_own, ID_nei);
                            scalar maxID = Foam::max(ID_own, ID_nei);


                            scalar diff(maxID - minID);
                            volScalarField filtered = diff
                                                  * pos( vofID_ - maxID + 0.1)
                                                  * pos(-vofID_ + maxID + 0.1);
                            vofID_ -= filtered;
                        }
                    }
                }
            }
        }
        vofID_.correctBoundaryConditions();
    }


    // Count number of continua and renumber vofID
    {
        volScalarField tmp(vofID_);
        vofID_ = 0.0;
        scalar maxID = gMax(tmp);
        int iter = 1;
        while (maxID > 0.1)
        {
            volScalarField curVolID = pos(tmp - maxID + 0.1);
            vofID_ += iter*curVolID;
            tmp -= maxID*curVolID;
            maxID = gMax(tmp);
            iter++;
        }
    }


    // Calculate bubble volume, velocity, and position
    // create lists to store data
    label maxID = floor(gMax(vofID_) + 0.5);
    labelList noCells(maxID+1, label(0));
    scalarList cellVolume(maxID+1, scalar(0));
    scalarList volume(maxID+1, scalar(0));
    vectorList position(maxID+1, vector(0,0,0));
    vectorList velocity(maxID+1, vector(0,0,0));

    
    // Loop over all cells and store corresponding data
    forAll(mesh_.cells(), cellI)
    {
        if (vofID_[cellI] > 0.1)
        {
            // Get volume ID
            label volID = floor(vofID_[cellI] + 0.5);

            // Store data in corresponding list
        	    
	    noCells[volID] += 1;
            cellVolume[volID] += mesh_.V()[cellI];
            volume[volID] += (1-alpha_[cellI])*mesh_.V()[cellI];
            velocity[volID] += (1-alpha_[cellI])*mesh_.V()[cellI]*U_[cellI];
            position[volID] += (1-alpha_[cellI])*mesh_.V()[cellI]*mesh_.C()[cellI];

       }
    }

    // Account for parallel processing
    reduce( noCells, sumOp<labelList>() );
    reduce( cellVolume, sumOp<scalarList>() );
    reduce( volume, sumOp<scalarList>() );
    reduce( position, sumOp<vectorList>() );
    reduce( velocity, sumOp<vectorList>() );
    // Calculate and inject bubble    (Vallier,2013) --- 每个网格注入一个bubble;
    //
    for(int i=1; i<=maxID; i++)
    {
	
        // correct the threshold values Nel;
        if ((volume[i] > VSMALL && volume[i] < vlim_ && VoFIDactive_) || (noCells[i] < Nel_ && NelCriterion_))
        {
            // Calculate bubble Diameter
            const scalar d = cbrt(6.0*volume[i]/constant::mathematical::pi);
	    
            // Weighting of position and velocity
            vector pos = position[i] / volume[i];
            vector U = velocity[i] / volume[i];
            
	    // Inject bubble
            cloud_.injectbubble(pos, d, U);

	    // Remove the Eulerian cavity of cell
            forAll(mesh_.cells(), cellI)
	        {
	            if( vofID_[cellI] == i )
		        {
			        beta_[cellI] = alpha_[cellI];
                                alpha_[cellI] = 1;
                                
		        }

	        }
           if(Pstream::master())
	   {
 		Info<<"Remove small Euler bubble and inject Lagrange bubble successful"<<endl;
 	   }	    
        }
    }
    // lagrangian bubble to VOF interface
}

// ************************************************************************* //
