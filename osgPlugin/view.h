#include "stdafx.h"
#include <osg/Group>
#include <osgViewer/Viewer>
#include <osgViewer/ViewerEventHandlers>
#include <osgViewer/api/win32/GraphicsWindowWin32>
#include <osgDB/ReadFile>
#include <osgGA/OrbitManipulator>
#include <osgGA/StateSetManipulator>
#include <osgGA/CameraManipulator>
#include <osgGA/KeySwitchMatrixManipulator>
#include <osgGA/GUIEventHandler>
#include <osgUtil/Optimizer>
#include <iostream>
#include <sstream>
#include <osgDB/ReadFile>
#include <osgViewer/Viewer>
#include <osg/Geode>
#include <osg/Depth>
#include <osg/CameraNode>
#include <osgText/Text>
#include <osgGA/TrackballManipulator>
#include <osg/LineWidth>
#include <osg/Point>
#include <osg/ShapeDrawable>
#include <osg/MatrixTransform>
#include <osgDB/ReadFile>
#include <osgGA/CameraManipulator>
#include <osgViewer/ViewerEventHandlers>
#include <osgDB/WriteFile>
#include <winsock2.h>
#include <osg/ShapeDrawable>
#include <osg/Shape>
using namespace std;
osg::ref_ptr<osg::Group> g_grpMouse;
char worldCoordinate[255];
class CHUD_viewPoint: public osgGA::GUIEventHandler  

{

public:  

 

    /**���캯��*/

      CHUD_viewPoint(osgText::Text* updateText):

      m_text(updateText) {}

 

      ~CHUD_viewPoint(){}

      virtual bool handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa);

      void UpdateText(osgViewer::Viewer* viewer,const osgGA::GUIEventAdapter&);

 

      /**LABEL*/

      void setLabel(const std::string& name)

      {

         if ( m_text.get())

         {

             m_text->setText(name); 

         }

      }
	  //void setLabel(const std::string& name)
	  //{
		 // if (_updateText.get()) _updateText->setText(name);
	  //}
	  float mX;
	  float mY;
	  bool mEnableDragger;

protected: 

    osg::Vec2 m_vPosWindowMouse;//��굥�����Ĵ�������

 	void pick(const osgGA::GUIEventAdapter& ea, osgGA::GUIActionAdapter& aa);

    osg::ref_ptr<osgText::Text>  m_text;//�ӵ���Ϣ���ᶯ̬�ı�

};
bool CHUD_viewPoint::handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa)

{

    switch(ea.getEventType())

    {

    case( osgGA::GUIEventAdapter::PUSH):

    {

       m_vPosWindowMouse.set( ea.getX(), ea.getY());//��굥�����Ĵ�������

      

       osgViewer::Viewer* viewer = dynamic_cast< osgViewer::Viewer*>( &aa);

      

       if (viewer&&!mEnableDragger)

       {

          // UpdateText( viewer, ea);//����������

           //�����

           osg::ref_ptr<osg::Camera> cameraMaster = viewer->getCamera();
           osg::Matrix mvpw = cameraMaster->getViewMatrix() * cameraMaster->getProjectionMatrix();

           if ( cameraMaster->getViewport()) mvpw.postMult( cameraMaster->getViewport()->computeWindowMatrix());

           osg::Matrix _inverseMVPW;

           _inverseMVPW.invert( mvpw);

           osg::Vec3d nearPoint = osg::Vec3d( ea.getX(), ea.getY(), 0.0)* _inverseMVPW;//͸��ͶӰ��Znearƽ��Ľ���

           osg::Vec3d farPoint = osg::Vec3d( ea.getX(), ea.getY(), 1.0)* _inverseMVPW;//͸��ͶӰ��Zfarƽ��Ľ���

           osg::Vec3 vPosEye, vCenter, vUp;

           cameraMaster->getViewMatrixAsLookAt( vPosEye, vCenter, vUp);//��ȡ�ӵ���Ϣ

           osg::Matrix _inverseMV;

           _inverseMV.invert( cameraMaster->getViewMatrix());

           osg::Vec3 ptEye= osg::Vec3(  0, 0, 0) * _inverseMV;//��ȡ�ӵ�����

           osg::Vec3d deltaEye= ptEye- vPosEye;

           if ( deltaEye.length()< 1e-8)

           {

              cout<< "yes,eye\n";

           }

           else

           {

              cout<< "no,eye\n";

           }

           osg::Vec3d dir1= farPoint- nearPoint;

           dir1.normalize();

           osg::Vec3d dir2= farPoint- vPosEye;

           dir2.normalize();

           osg::Vec3d delta= dir1- dir2;

           //���ӵ㡢Znearƽ��Ľ��㡢Zfarƽ��Ľ����Ƿ���ͬһֱ���ϡ�����֤��ȷ����ͬһֱ����

           if ( delta.length()< 1e-8)

           {

              cout<< "yes,line\n";

           }

           else

           {

              cout<< "no,line\n";

           }

           osg::Geode* geode= new osg::Geode();

           osg::Geometry* pyramidGeometry = new osg::Geometry();

           geode->addDrawable( pyramidGeometry);

           osg::Vec3Array* pyramidVertices = new osg::Vec3Array;

           pyramidVertices->push_back( nearPoint);

           pyramidVertices->push_back( farPoint);

           pyramidGeometry->setVertexArray( pyramidVertices );

           //��ɫ

           osg::Vec4Array* colors = new osg::Vec4Array;

           colors->push_back( osg::Vec4(  1.0f, 0.0f, 0.0f, 1.0f) );//��ɫ

           pyramidGeometry->setColorArray( colors);

           pyramidGeometry->setColorBinding( osg::Geometry::BIND_OVERALL);

          

           //����ʾ͸��ͶӰ��Znearƽ��Ľ���

          pyramidGeometry->addPrimitiveSet(  new osg::DrawArrays( osg::PrimitiveSet::POINTS, 0, 1/*3*/));

           //���߱�ʾ��������ߣ������ΪZnearƽ�潻�㣬�յ�ΪZfarƽ�潻�㡣

           pyramidGeometry->addPrimitiveSet(  new osg::DrawArrays( osg::PrimitiveSet::LINES, 0, 2));/**/


           osg::ref_ptr <osg::Point> ptSize = new osg::Point;

           ptSize->setSize( 0.5) ;      

           geode->getOrCreateStateSet()->setAttributeAndModes( ptSize.get (),osg::StateAttribute::ON);          

          

           /*��ֻ��һ����ʱ����Χ��뾶Ϊ�����Կ��ܿ���������㣬����Ҫ�������ð�Χ���С���ɰѰ�Χ��뾶���㡣

           ���glider��cow��Сģ�ͣ��뾶ȡ.1���ԣ���fountain.osg��.1̫С��Ϊͳһ���ɴ�Щ����*/

           osg::Vec3d ptCnt= geode->getBound().center();

           double dRadius= geode->getBound().radius();

           //�������ð�Χ��İ뾶���ɵ���setInitialBound()��

           osg::BoundingSphere bs( ptCnt, 0.1);

           geode->setInitialBound( bs);

 
           g_grpMouse->removeChildren( 0, g_grpMouse->getNumChildren());

           g_grpMouse->addChild( geode);


           //��ȡ�Ӹ��ڵ㵽��ǰ�ڵ��·������

           osg::NodePathList parentNodePaths = geode->getParentalNodePaths();

           if ( !parentNodePaths.empty())
           {
              osg::Matrixd mt= computeWorldToLocal( parentNodePaths[ 0]);
           }

           osg::ref_ptr< osgUtil::LineSegmentIntersector > picker =new osgUtil::LineSegmentIntersector(

              nearPoint, farPoint);//�߶�(��ʵ����������)

           osgUtil::IntersectionVisitor iv( picker.get());

           g_grpMouse->getParent( 0)->getChild( 0)->/*asGroup()->getChild( 0)->*/accept( iv);//ģ����/**/

           if (picker->containsIntersections())

           {

              osg::Vec3 ptWorldIntersectPointFirst= picker->getFirstIntersection().getWorldIntersectPoint();

              cout<<"world coords vertex("<< ptWorldIntersectPointFirst.x()<<","

                  << ptWorldIntersectPointFirst.y()<< ","<< ptWorldIntersectPointFirst.z()<<")"<< std::endl;

			  std::string gdlist=""; 

			  std::ostringstream os;

			  os<<"World Coordinate:("<< ptWorldIntersectPointFirst.x()<<","

				  << ptWorldIntersectPointFirst.y()<< ","<< ptWorldIntersectPointFirst.z()<<")";//����

			  gdlist = os.str(); 
			//  strCoordinate = gdlist.c_str();;
			  setLabel(gdlist); 

			  std::ostringstream os2;
			  os2<<ptWorldIntersectPointFirst.x()<<","

				  << ptWorldIntersectPointFirst.y()<< ","<< ptWorldIntersectPointFirst.z();//����
			  gdlist = os2.str(); 
			  strcpy(worldCoordinate,gdlist.c_str());
              //�����˵�              

              double dPointRadius= 0.5f;

              osg::ShapeDrawable* pShd= new osg::ShapeDrawable(

                  new osg::Sphere( ptWorldIntersectPointFirst, dPointRadius));

              pShd->setColor( osg::Vec4( 0, 1, 0, 1));

              geode->addDrawable( pShd);

           }
       }
	   else
	   {
		   mX = ea.getX();
		   mY = ea.getY();
	   }

      	break;
		//return true;

    }
	case osgGA::GUIEventAdapter::RELEASE:
		{
			if (ea.getX() == mX && ea.getY() == mY)
			{
				pick(ea, aa);
			}
		}
		break;

	case osgGA::GUIEventAdapter::KEYDOWN:
		{
			if (ea.getKey() == 'd')
			{
				mEnableDragger = !mEnableDragger;
			}
			if (ea.getKey() == 'f'||ea.getKey() == 'F')
			{
				
			}
			if (ea.getKey() == 0x20)
			{
				aa.requestRedraw();
				aa.requestContinuousUpdate(false);
			}
		}
		break;
    default:     

       return false; 
    }

}
void CHUD_viewPoint::pick(const osgGA::GUIEventAdapter &ea, osgGA::GUIActionAdapter &aa)
{
	osgViewer::View* view = dynamic_cast<osgViewer::View*> (&aa); 

	osgUtil::LineSegmentIntersector::Intersections hits;
	if (view->computeIntersections(ea.getX(), ea.getY(), hits))
	{
		osgUtil::LineSegmentIntersector::Intersection intersection = *hits.begin();
		osg::NodePath& nodePath = intersection.nodePath;
		int nNodeSize = static_cast<int> (nodePath.size());

		if (nNodeSize > 0)
		{
			osg::Node* node = nodePath[nNodeSize - 1];
			osg::Node* grandParent = node->getParent(0)->getParent(0);

			// This method maybe not right?
			ModelShape* shape = dynamic_cast<ModelShape*> (grandParent);
			if (shape)
			{
				mEnableDragger ? shape->EnableDragger() : shape->DisableDragger();
			}

		}
	}

}
void CHUD_viewPoint::UpdateText(osgViewer::Viewer* viewer,const osgGA::GUIEventAdapter&)

{

    std::string gdlist=""; 

 

    std::ostringstream os;

    os<<"MousePos(X: "<< m_vPosWindowMouse.x()<<",Y: "<< m_vPosWindowMouse.y()<<")";//����

   

    gdlist = os.str(); 

 

    setLabel(gdlist); 

}
 osg::Node* createHUD_viewPoint( osgText::Text* text)

{

    //��������

    std::string font("fonts/arial.TTF");//�˴����õ��Ǻ������� "fonts/STCAIYUN.TTF"

    text->setFont( font); 

    //����������ʾ��λ��(����Ϊ(0,0),X�����ң�Y������)

    osg::Vec3 position( 100.0f, 10.0f,0.0f);

    text->setPosition(position);   

    text->setColor( osg::Vec4( 1, 1, 0, 1));

    text->setText(L"");//������ʾ������

    text->setCharacterSize(15);

    text->setDataVariance(osg::Object::DYNAMIC);//һ��Ҫ��������Ϊ��̬���������Ῠס���������������osgcatch��

 

    //������ڵ�

    osg::Geode* geode = new osg::Geode();

    geode->addDrawable( text );//������Text����drawable���뵽Geode�ڵ���

    //����״̬

    osg::StateSet* stateset = geode->getOrCreateStateSet(); 

    stateset->setMode(GL_LIGHTING,osg::StateAttribute::OFF);//�رյƹ�

    stateset->setMode(GL_DEPTH_TEST,osg::StateAttribute::OFF);//�ر���Ȳ���

    //��GL_BLEND���ģʽ���Ա�֤Alpha������ȷ��

    stateset->setMode(GL_BLEND,osg::StateAttribute::ON);

 

    //���

    osg::Camera* camera = new osg::Camera;

    //����͸�Ӿ���

    camera->setProjectionMatrix(osg::Matrix::ortho2D(0,600,0,600));//����ͶӰ   

    //���þ��Բο�����ϵ��ȷ����ͼ���󲻻ᱻ�ϼ��ڵ�ı任����Ӱ��

    camera->setReferenceFrame(osg::Transform::ABSOLUTE_RF);

    //��ͼ����ΪĬ�ϵ�

    camera->setViewMatrix(osg::Matrix::identity());

 

    //���ñ���Ϊ͸��������Ļ���������ClearColor 

    camera->setClearMask(GL_DEPTH_BUFFER_BIT);

    camera->setAllowEventFocus( false);//����Ӧ�¼���ʼ�յò�������

 

    //������Ⱦ˳�򣬱����������Ⱦ

    camera->setRenderOrder(osg::CameraNode::POST_RENDER); 

 

    camera->addChild(geode);//��Ҫ��ʾ��Geode�ڵ���뵽���

 

    return camera; 

};