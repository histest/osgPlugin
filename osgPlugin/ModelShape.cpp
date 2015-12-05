#include <stdafx.h>
#include "ModelShape.h"

ModelShape::ModelShape(osg::Node* shape)
	: mShape(shape)
	, mDragger(new osgManipulator::TrackballDragger())
	,mDragger2(new osgManipulator::TranslateAxisDragger())
	, mSelection(new osgManipulator::Selection())
{
	float scale = shape->getBound().radius() * 0.2;
	mDragger->setMatrix(osg::Matrix::scale(scale, scale, scale) * 
		osg::Matrix::translate(shape->getBound().center()));

	mDragger->setupDefaultGeometry();

	mSelection->addChild(shape);

	scale = shape->getBound().radius() * 1.5;
	mDragger2->setMatrix(osg::Matrix::scale(scale, scale, scale) * 
		osg::Matrix::translate(shape->getBound().center()));

	mDragger2->setupDefaultGeometry();


	addChild(mSelection);
}


ModelShape::~ModelShape(void)
{
}


void ModelShape::EnableDragger()
{
	addChild(mDragger);
	addChild(mDragger2);
	mDragger->addTransformUpdating(mSelection);
	mDragger->setHandleEvents(true);

	mDragger2->addTransformUpdating(mSelection);
	mDragger2->setHandleEvents(true);
}


void ModelShape::DisableDragger()
{
	removeChild(mDragger);

	mDragger->removeTransformUpdating(mSelection);
	mDragger->setHandleEvents(false);

	removeChild(mDragger2);

	mDragger2->removeTransformUpdating(mSelection);
	mDragger2->setHandleEvents(false);
}