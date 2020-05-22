//
//  Scenes.h
//  MultiTracker
//
//  Created by Fang Da on 15/1/26.
//
//

#ifndef __MultiTracker__Scenes__
#define __MultiTracker__Scenes__

#include <iostream>
#include "VS3D.h"

class Scenes
{
public:
    // scene-specific initialization
    static VS3D * sceneSphere              ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneTet                 ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneCube                ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneSheet               ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneBarrel              ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneDoubleBubble        ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneTwoBubbles          ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneTripleJunction      ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneFoamInit            ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneFoam                ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneQuadJunction        ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneConstrainedSphere   ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneBubbleWand          ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneTwoRingsPinching    ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * scenePullingFoam         ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * scenePeanutBubble        ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneStraw               ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneCarousel            ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneOctahedron          ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneBubbleLattice       ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneMergedBubbleLattice ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneFlyingBubbles       ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneBubbleLine          ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneBlowingBubble          ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneStraws          ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * sceneNewFoam          ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static VS3D * scene2DBubbleLattice          ( Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    
    // scene-specific time evolution
    static void stepSphere              ( double dt, Sim * sim, VS3D * vs);
    static void stepTet                 ( double dt, Sim * sim, VS3D * vs);
    static void stepCube                ( double dt, Sim * sim, VS3D * vs);
    static void stepSheet               ( double dt, Sim * sim, VS3D * vs);
    static void stepBarrel              ( double dt, Sim * sim, VS3D * vs);
    static void stepDoubleBubble        ( double dt, Sim * sim, VS3D * vs);
    static void stepTwoBubbles          ( double dt, Sim * sim, VS3D * vs);
    static void stepTripleJunction      ( double dt, Sim * sim, VS3D * vs);
    static void stepFoamInit            ( double dt, Sim * sim, VS3D * vs);
    static void stepFoam                ( double dt, Sim * sim, VS3D * vs);
    static void stepQuadJunction        ( double dt, Sim * sim, VS3D * vs);
    static void stepConstrainedSphere   ( double dt, Sim * sim, VS3D * vs);
    static void stepBubbleWand          ( double dt, Sim * sim, VS3D * vs);
    static void stepTwoRingsPinching    ( double dt, Sim * sim, VS3D * vs);
    static void stepPullingFoam         ( double dt, Sim * sim, VS3D * vs);
    static void stepPeanutBubble        ( double dt, Sim * sim, VS3D * vs);
    static void stepStraw               ( double dt, Sim * sim, VS3D * vs);
    static void stepCarousel            ( double dt, Sim * sim, VS3D * vs);
    static void stepOctahedron          ( double dt, Sim * sim, VS3D * vs);
    static void stepBubbleLattice       ( double dt, Sim * sim, VS3D * vs);
    static void stepMergedBubbleLattice ( double dt, Sim * sim, VS3D * vs);
    static void stepFlyingBubbles       ( double dt, Sim * sim, VS3D * vs);
    static void stepBubbleLine          ( double dt, Sim * sim, VS3D * vs);
    static void stepBlowingBubble          ( double dt, Sim * sim, VS3D * vs);
    static void stepStraws          ( double dt, Sim * sim, VS3D * vs);
    static void stepNewFoam          ( double dt, Sim * sim, VS3D * vs);
    static void step2DBubbleLattice          ( double dt, Sim * sim, VS3D * vs);

};

#endif /* defined(__MultiTracker__Scenes__) */
