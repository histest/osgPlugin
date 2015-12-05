// ConsoleApplication1.cpp : 定义控制台应用程序的入口点。
//


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
 

#pragma comment( lib, "osgd.lib"); //.在Debug版本下的库名都加d,如"osgd.lib"

#pragma comment( lib, "osgDBd.lib")

#pragma comment( lib, "osgViewerd.lib");

#pragma comment( lib, "osgTextd.lib");

#pragma comment( lib, "osgGAd.lib");

#pragma comment( lib, "osgUtild.lib");

 

osg::ref_ptr<osg::Group> g_grpMouse;

using namespace std;

 char worldCoordinate[255];
 

//事件类

class CHUD_viewPoint: public osgGA::GUIEventHandler  

{

public:  

 

    /**构造函数*/

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
 

protected: 

    osg::Vec2 m_vPosWindowMouse;//鼠标单击处的窗口坐标

 

    osg::ref_ptr<osgText::Text>  m_text;//视点信息，会动态改变

};

 bool CHUD_viewPoint::handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa)

{

    switch(ea.getEventType())

    {

    //case(osgGA::GUIEventAdapter::FRAME):

    //  {

    //     osgViewer::Viewer* viewer = dynamic_cast<osgViewer::Viewer*>(&aa);

 

    //     if (viewer)

    //     {

    //         osg::Vec3 vCenter, vUp;

    //         viewer->getCamera()->getViewMatrixAsLookAt( m_vPosEye, vCenter, vUp);//获取视点信息

    //         UpdateText( viewer, ea);//更新文字信息       

 

    //     }

 

    //     return true; 

    //  }    

    case( osgGA::GUIEventAdapter::PUSH):

    {

       m_vPosWindowMouse.set( ea.getX(), ea.getY());//鼠标单击处的窗口坐标

      

       osgViewer::Viewer* viewer = dynamic_cast< osgViewer::Viewer*>( &aa);

      

       if (viewer)

       {

          // UpdateText( viewer, ea);//更新文字信息

 

           //主相机

           osg::ref_ptr<osg::Camera> cameraMaster = viewer->getCamera();

 

           osg::Matrix mvpw = cameraMaster->getViewMatrix() * cameraMaster->getProjectionMatrix();

           if ( cameraMaster->getViewport()) mvpw.postMult( cameraMaster->getViewport()->computeWindowMatrix());

 

           osg::Matrix _inverseMVPW;

           _inverseMVPW.invert( mvpw);

 

           osg::Vec3d nearPoint = osg::Vec3d( ea.getX(), ea.getY(), 0.0)* _inverseMVPW;//透视投影中Znear平面的交点

           osg::Vec3d farPoint = osg::Vec3d( ea.getX(), ea.getY(), 1.0)* _inverseMVPW;//透视投影中Zfar平面的交点

 

           osg::Vec3 vPosEye, vCenter, vUp;

           cameraMaster->getViewMatrixAsLookAt( vPosEye, vCenter, vUp);//获取视点信息

 

 

           osg::Matrix _inverseMV;

           _inverseMV.invert( cameraMaster->getViewMatrix());

           osg::Vec3 ptEye= osg::Vec3(  0, 0, 0) * _inverseMV;//获取视点坐标

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

 

           //看视点、Znear平面的交点、Zfar平面的交点是否在同一直线上。经验证，确定在同一直线上

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

 

           //颜色

           osg::Vec4Array* colors = new osg::Vec4Array;

           colors->push_back( osg::Vec4(  1.0f, 0.0f, 0.0f, 1.0f) );//红色

           pyramidGeometry->setColorArray( colors);

           pyramidGeometry->setColorBinding( osg::Geometry::BIND_OVERALL);

          

           //红点表示透视投影中Znear平面的交点

          pyramidGeometry->addPrimitiveSet(  new osg::DrawArrays( osg::PrimitiveSet::POINTS, 0, 1/*3*/));

           //红线表示鼠标点击的线，其起点为Znear平面交点，终点为Zfar平面交点。

           pyramidGeometry->addPrimitiveSet(  new osg::DrawArrays( osg::PrimitiveSet::LINES, 0, 2));/**/

 

           ////设置线宽

           //osg::ref_ptr <osg::LineWidth> LineSize = new osg::LineWidth;

           //LineSize ->setWidth( 12.0) ;      

           //geode->getOrCreateStateSet()->setAttributeAndModes( LineSize.get (),osg::StateAttribute::ON);

 

           //设置点大小

           osg::ref_ptr <osg::Point> ptSize = new osg::Point;

           ptSize->setSize( 0.5) ;      

           geode->getOrCreateStateSet()->setAttributeAndModes( ptSize.get (),osg::StateAttribute::ON);          

          

           /*当只有一个点时，包围球半径为，所以可能看不到这个点，故需要重新设置包围球大小，可把包围球半径设大点。

           如对glider、cow等小模型，半径取.1可以，对fountain.osg则.1太小。为统一，可大些，如*/

           osg::Vec3d ptCnt= geode->getBound().center();

           double dRadius= geode->getBound().radius();

           //重新设置包围球的半径（可调用setInitialBound()）

           osg::BoundingSphere bs( ptCnt, 10);

           geode->setInitialBound( bs);

 

           g_grpMouse->removeChildren( 0, g_grpMouse->getNumChildren());

           g_grpMouse->addChild( geode);

 

           //获取从根节点到当前节点的路径向量

           osg::NodePathList parentNodePaths = geode->getParentalNodePaths();

           if ( !parentNodePaths.empty())

           {

              osg::Matrixd mt= computeWorldToLocal( parentNodePaths[ 0]);

 

           }

 

           ////求交

           //osg::ref_ptr< osgUtil::LineSegmentIntersector > picker = new osgUtil::LineSegmentIntersector(

           //  osgUtil::Intersector::WINDOW, ea.getX(), ea.getY());

           //osgUtil::IntersectionVisitor iv( picker.get());

           ////g_grpMouse->getParent( 0)->getChild( 0)->accept( iv);//模型求交

           //cameraMaster->accept( iv);//模型求交(从相机往下遍历)

 

           //求交

           osg::ref_ptr< osgUtil::LineSegmentIntersector > picker =new osgUtil::LineSegmentIntersector(

              nearPoint, farPoint);//线段(真实的世界坐标)

           osgUtil::IntersectionVisitor iv( picker.get());

           //g_grpMouse->getParent( 0)->getChild( 0)->accept( iv);//模型求交/**/

           g_grpMouse->getParent( 0)->getChild( 0)->/*asGroup()->getChild( 0)->*/accept( iv);//模型求交/**/

            //               根节点     cow的MT节点

           //cameraMaster->accept( iv);//模型求交(从相机往下遍历)

 

 

          

           //if (picker->containsIntersections())

           //{

           //  double dLen2Shortest= DBL_MAX, dLenTmp;

           //  osgUtil::LineSegmentIntersector::Intersections::iterator hitrShortest;

           //  osgUtil::LineSegmentIntersector::Intersections intersections= picker->getIntersections();

 

           //  for( osgUtil::LineSegmentIntersector::Intersections::iterator hitr = intersections.begin();

           //     hitr != intersections.end();

           //     ++hitr)

           //  {

           //     //求离视点最近的点，即鼠标选择的最前面的点

           //     dLenTmp= ( ptEye.x()- hitr->getWorldIntersectPoint().x())*

           //         ( ptEye.x()- hitr->getWorldIntersectPoint().x())+

           //         ( ptEye.y()- hitr->getWorldIntersectPoint().y())*

           //         ( ptEye.y()- hitr->getWorldIntersectPoint().y())+

           //         ( ptEye.z()- hitr->getWorldIntersectPoint().z())*

           //         ( ptEye.z()- hitr->getWorldIntersectPoint().z());

           //     if ( dLen2Shortest> dLenTmp)

           //     {

           //         dLen2Shortest= dLenTmp;

           //         hitrShortest= hitr;                   

           //     }

           //  }

           // 

           //  //输出

           //  if ( dLen2Shortest != DBL_MAX)

           //  {

           //     cout<<"world coords vertex("<< hitrShortest->getWorldIntersectPoint().x()<<","

           //         << hitrShortest->getWorldIntersectPoint().y()<<","

           //         << hitrShortest->getWorldIntersectPoint().z()<<")"<<std::endl;

 

           //     //高亮此点              

           //     double dPointRadius= 15.0f;

           //     osg::ShapeDrawable* pShd= new osg::ShapeDrawable(

           //         new osg::Sphere( hitrShortest->getWorldIntersectPoint(), dPointRadius));//绘制交点的球

           //     pShd->setColor( osg::Vec4( 0, 1, 0, 1));

           //     geode->addDrawable( pShd);

           //    

           //  }

 

 

           //}

 

           if (picker->containsIntersections())

           {

              osg::Vec3 ptWorldIntersectPointFirst= picker->getFirstIntersection().getWorldIntersectPoint();

              cout<<"world coords vertex("<< ptWorldIntersectPointFirst.x()<<","

                  << ptWorldIntersectPointFirst.y()<< ","<< ptWorldIntersectPointFirst.z()<<")"<< std::endl;

			  std::string gdlist=""; 



			  std::ostringstream os;

			  os<<"World Coordinate:("<< ptWorldIntersectPointFirst.x()<<","

				  << ptWorldIntersectPointFirst.y()<< ","<< ptWorldIntersectPointFirst.z()<<")";//坐标



			  gdlist = os.str(); 
			//  strCoordinate = gdlist.c_str();;
			  strcpy(worldCoordinate,gdlist.c_str());

			  setLabel(gdlist); 

              //高亮此点              

              double dPointRadius= 1.0f;

              osg::ShapeDrawable* pShd= new osg::ShapeDrawable(

                  new osg::Sphere( ptWorldIntersectPointFirst, dPointRadius));

              pShd->setColor( osg::Vec4( 0, 1, 0, 1));

              geode->addDrawable( pShd);

 

           }

          

 

       }

       return true;

    }

    default:     

       return false; 

 

    }

}

void CHUD_viewPoint::UpdateText(osgViewer::Viewer* viewer,const osgGA::GUIEventAdapter&)

{

    std::string gdlist=""; 

 

    std::ostringstream os;

    os<<"MousePos(X: "<< m_vPosWindowMouse.x()<<",Y: "<< m_vPosWindowMouse.y()<<")";//坐标

   

    gdlist = os.str(); 

 

    setLabel(gdlist); 

}

 osg::Node* createHUD_viewPoint( osgText::Text* text)

{

    //设置字体

    std::string font("fonts/arial.TTF");//此处设置的是汉字字体 "fonts/STCAIYUN.TTF"

    text->setFont( font); 

    //设置文字显示的位置(左下为(0,0),X正向朝右，Y正向朝上)

    osg::Vec3 position( 100.0f, 10.0f,0.0f);

    text->setPosition(position);   

    text->setColor( osg::Vec4( 1, 1, 0, 1));

    text->setText(L"");//设置显示的文字

    text->setCharacterSize(15);

    text->setDataVariance(osg::Object::DYNAMIC);//一定要设置字体为动态，否则程序会卡住，死在那里。（参照osgcatch）

 

    //几何体节点

    osg::Geode* geode = new osg::Geode();

    geode->addDrawable( text );//将文字Text作这drawable加入到Geode节点中

    //设置状态

    osg::StateSet* stateset = geode->getOrCreateStateSet(); 

    stateset->setMode(GL_LIGHTING,osg::StateAttribute::OFF);//关闭灯光

    stateset->setMode(GL_DEPTH_TEST,osg::StateAttribute::OFF);//关闭深度测试

    //打开GL_BLEND混合模式（以保证Alpha纹理正确）

    stateset->setMode(GL_BLEND,osg::StateAttribute::ON);

 

    //相机

    osg::Camera* camera = new osg::Camera;

    //设置透视矩阵

    camera->setProjectionMatrix(osg::Matrix::ortho2D(0,600,0,600));//正交投影   

    //设置绝对参考坐标系，确保视图矩阵不会被上级节点的变换矩阵影响

    camera->setReferenceFrame(osg::Transform::ABSOLUTE_RF);

    //视图矩阵为默认的

    camera->setViewMatrix(osg::Matrix::identity());

 

    //设置背景为透明，否则的话可以设置ClearColor 

    camera->setClearMask(GL_DEPTH_BUFFER_BIT);

    camera->setAllowEventFocus( false);//不响应事件，始终得不到焦点

 

    //设置渲染顺序，必须在最后渲染

    camera->setRenderOrder(osg::CameraNode::POST_RENDER); 

 

    camera->addChild(geode);//将要显示的Geode节点加入到相机

 

    return camera; 

};


#include <osgDB/ReadFile>
#include <osgGA/CameraManipulator>
#include <osgViewer/ViewerEventHandlers>
#include <osgDB/WriteFile>
#include <winsock2.h>
#include <osg/ShapeDrawable>
#include <osg/Shape>
#include "DataDistribute.h"
#include "Satellite.h"
#include "RTI_NET.h"
#include "OSGBUSChannelNo.h"
#include "Manipulator.h"
#include "Hand.h"
#include "CreateDevice.h"
#include <osg/Camera>
#include "FederalNo.h"
#include <conio.h>

osg::ref_ptr<osgViewer::Viewer> viewer; 
osg::ref_ptr<osg::Group> root;
osg::ref_ptr< osg::Node> model;
osg::ref_ptr< osg::Node> model2;
double TransformMatrix[3][3];
MainBody HumanStones(root.get(),viewer.get());		
 void CollectData(int Scale)
{
//	clsDevices[0].pMatrixTransform->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(45.0), 1, 0, 0) * osg::Matrix::scale(3, 3, 3) * osg::Matrix::translate(4, 0, -2));
	
	double scaled=Scale;
	HumanStones.clsActiveDevices[0].bRotateByQuaternion=false;
	HumanStones.clsActiveDevices[0].vEulerAttitude._v[0] = 10+scaled/10;
	HumanStones.clsActiveDevices[0].vEulerAttitude._v[1] = 10;
	HumanStones.clsActiveDevices[0].vEulerAttitude._v[2] = 10;
	
	HumanStones.clsActiveDevices[0].vPosition._v[0] = 100;
	HumanStones.clsActiveDevices[0].vPosition._v[1] = -40;//-1500*sin(scaled)
	HumanStones.clsActiveDevices[0].vPosition._v[2] = 300;
	HumanStones.clsActiveDevices[0].vScale.x()=2;
	HumanStones.clsActiveDevices[0].vScale.y()=2;
	HumanStones.clsActiveDevices[0].vScale.z()=2;

	HumanStones.clsActiveDevices[1].bRotateByQuaternion=false;
	HumanStones.clsActiveDevices[1].vEulerAttitude._v[0] = 10+scaled/10;;
	HumanStones.clsActiveDevices[1].vEulerAttitude._v[1] = 10;
	HumanStones.clsActiveDevices[1].vEulerAttitude._v[2] = 10;
	HumanStones.clsActiveDevices[1].vPosition._v[0] = 100;
	HumanStones.clsActiveDevices[1].vPosition._v[1] = 40;//+1500*sin(scaled)
	HumanStones.clsActiveDevices[1].vPosition._v[2] = 300;
	HumanStones.clsActiveDevices[1].vScale.x()=2;
	HumanStones.clsActiveDevices[1].vScale.y()=2;
	HumanStones.clsActiveDevices[1].vScale.z()=2;

 }
 /*char   buff[255]; */   
 const char* GetWorldIntersectPoint()
 {
	/* strcpy(buff,"100");*/
	 return worldCoordinate;
 }
#include "UserMechanics.h"
void  GetTransformMatrix(const double R[3],const double V[3],const double R2[3],const double V2[3],double Rbe[3][3])
 {
	 double Rbi[3][3] = {0.0},Rei[3][3] = {0.0};
	 GetRdi(R,V,Rbi);
	 GetRdi(R2,V2,Rei);
	 MatrixProductTranspose(3,3,3,&Rbi[0][0],&Rei[0][0],&TransformMatrix[0][0]);
	 for (int i=0;i<3;i++)
	 {
		 for (int j=0;j<3;j++)
		 {
			 Rbe[i][j]=TransformMatrix[i][j];
		 }
	 }
 }
int osgRunCow( HWND m_hWnd,char * Path,char*Path2)

{
	RECT rect;
	::GetWindowRect(m_hWnd, &rect);
	root= new osg::Group; 
	viewer=new osgViewer::Viewer;	

	HumanStones.clsActiveDevices[0].pNodeModel = osgDB::readNodeFile(Path);
	HumanStones.clsDevices[0].pNodeModel->getOrCreateStateSet()->setMode(GL_NORMALIZE,osg::StateAttribute::ON);
	HumanStones.clsActiveDevices[0].pMatrixTransform->addChild(HumanStones.clsActiveDevices[0].pNodeModel);
	HumanStones.clsActiveDevices[0].pMatrixTransform->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(45.0), 1, 0, 0) * osg::Matrix::scale(1, 1, 1) * osg::Matrix::translate(4, 0, -2));
	HumanStones.pMatrixTransform->addChild(HumanStones.clsActiveDevices[0].pMatrixTransform);
	HumanStones.clsActiveDevices[0].pMatrixTransform->setUpdateCallback(&(HumanStones.clsActiveDevices[0]));
	//root->addChild(HumanStones.pMatrixTransform);


	HumanStones.clsActiveDevices[1].pNodeModel = osgDB::readNodeFile(Path2);
	HumanStones.clsDevices[1].pNodeModel->getOrCreateStateSet()->setMode(GL_NORMALIZE,osg::StateAttribute::ON);
	HumanStones.clsActiveDevices[1].pMatrixTransform->addChild(HumanStones.clsActiveDevices[1].pNodeModel);
	HumanStones.clsActiveDevices[1].pMatrixTransform->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(145.0), 1, 0, 0) * osg::Matrix::scale(2, 2, 2) * osg::Matrix::translate(4, 0, -2));
	HumanStones.pMatrixTransform->addChild(HumanStones.clsActiveDevices[1].pMatrixTransform);
	HumanStones.clsActiveDevices[1].pMatrixTransform->setUpdateCallback(&(HumanStones.clsActiveDevices[1]));
	root->addChild(HumanStones.pMatrixTransform);



	model = osgDB::readNodeFile(Path);// glider nathan
	
   
	model2 = osgDB::readNodeFile(Path2);
	if (model!=NULL&&model2!=NULL)
	{
		/*root->addChild( model.get());

		root->addChild( model2.get());
		osg::ref_ptr<osg::MatrixTransform> trans = new osg::MatrixTransform;
		trans->setMatrix(osg::Matrix::translate(2000, 0, 0));
		trans->addChild(model2.get());

		osg::ref_ptr<osg::MatrixTransform> rot = new osg::MatrixTransform;
		rot->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(45.0), 1, 0, 0) * osg::Matrix::scale(1, 1, 1) * osg::Matrix::translate(4, 0, -2));
		rot->addChild(model2.get());
		root->addChild(trans.get());

		osg::ref_ptr<osg::PositionAttitudeTransform> pat = 
		new osg::PositionAttitudeTransform();
		pat->addChild(model);

		pat->setUpdateCallback(new RotateCallBack() );
		root ->addChild (pat) ;


		osg::MatrixTransform* transform = new osg::MatrixTransform();
		transform->setUpdateCallback(new MyTransformCallback(10.0f));
		transform->addChild(model2);

		root ->addChild (transform) ;*/
	}
  

    osgText::Text* text = new osgText::Text; 

    root->addChild( createHUD_viewPoint( text));//加入HUD文字

   

    osg::ref_ptr< CHUD_viewPoint> pHUD= new CHUD_viewPoint( text);

    viewer->addEventHandler( pHUD.get());

 

 
	// 初始化图形环境GraphicsContext Traits

	osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits;
	// Init the Windata Variable that holds the handle for the Window to display OSG in.
	osgViewer::GraphicsWindowWin32::WindowData* windata = new osgViewer::GraphicsWindowWin32::WindowData(m_hWnd);
	// 设置traits 参数
	traits->x = 0;
	traits->y = 0;
	traits->width = rect.right - rect.left;
	traits->height = rect.bottom - rect.top;
	traits->windowDecoration = false;
	traits->doubleBuffer = true;
	traits->sharedContext = 0;
	traits->setInheritedWindowPixelFormat = true;
	traits->inheritedWindowData = windata;
	//// 创建图形环境Graphics Context
	 osg::ref_ptr<osg::GraphicsContext> gc = osg::GraphicsContext::createGraphicsContext(traits.get());
	{
		// 根据分辨率确定合适的投影来保证显示的图形不变形
		/*double fovy =0.0; 
		double Ratio =0.0;
		double zNear =0.0;
		double zFar= 0.0;
		viewer.getCamera()->getProjectionMatrixAsPerspective(fovy,Ratio,zNear,zFar);
		double newAspectRatio =double(traits->width)/double(traits->height);
		double aspectRatioChange = newAspectRatio/aspectRatio;
		if (aspectRatioChange != 1.0)
		{
		viewer.getCamera()->getProjectionMatrix() *= osg::Matrix::scale(1.0/aspectRatioChange,1.0,1.0);
		}*/
	}


    gc->setClearColor( osg::Vec4f( 0.0f, 1.0f, 0.0f, 1.0f)); //设置整个windows窗口颜色

    gc->setClearMask( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   

    //主相机

    osg::ref_ptr<osg::Camera> cameraMaster = viewer->getCamera();

    cameraMaster->setGraphicsContext(gc.get());//设置GraphicsContext设备上下文

    //相机视口设置

    cameraMaster->setViewport(new osg::Viewport( 0, 0, traits->width, traits->height));/**/
   

    g_grpMouse= new osg::Group();

    //设置状态

    osg::StateSet* stateset = g_grpMouse->getOrCreateStateSet(); 

    stateset->setMode( GL_LIGHTING,osg::StateAttribute::OFF);//关闭灯光

    //stateset->setMode( GL_DEPTH_TEST,osg::StateAttribute::OFF);//关闭深度测试

    root->addChild( g_grpMouse.get());

 

    //viewer.setUpViewInWindow( 0, 0, 200, 200);//设置窗口大小


	DWORD dwThreadID;
	HANDLE hThread;

#ifdef CONNECT_ARCHIVE
//	hThread = CreateThread(NULL,0,(LPTHREAD_START_ROUTINE)CollectData,0,0,&dwThreadID);  /*USE CHANNEL ulDBArchiveChannelNoforArchiveAgent */
#endif



    viewer->setSceneData( root.get());

    viewer->realize();

    viewer->run() ;  
	
   
    return 0;

}