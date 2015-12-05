// osgPlugin.cpp : 定义 DLL 应用程序的导出函数。
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
#include <osgDB/ReadFile>
#include <osgGA/CameraManipulator>
#include <osgViewer/ViewerEventHandlers>
#include <osgDB/WriteFile>
#include <osg/ShapeDrawable>
#include <osg/Shape>
#include <osg/Camera>
#include <conio.h>
#include "UserMechanics.h"
#include "ModelShape.h"
#include "PickHandler.h"
#include "UnitConstant.h"
#include <osg/PositionAttitudeTransform>
#include "view.h"

osg::ref_ptr<osgViewer::Viewer> viewer; 
osg::ref_ptr<osg::Group> root;
osg::ref_ptr< osg::Node> model;
osg::ref_ptr< osg::Node> model2;
static double Rpt0[3][3]={0.0};//骨头上点在tracker系姿态
static double Rsp_right[3][3]={0.0};//模型上点在模型坐标系下姿态(固定值)
static double Rot0[3][3]={0.0};//
static double Rbo0[3][3]={0.0};//
static double Rpo0[3][3]={0.0};//模型上三点在osg下姿态

static double Rpt[3][3]={0.0};//骨头上点在tracker系姿态
static double Rsp_left[3][3]={0.0};//模型上点在模型坐标系下姿态(固定值)
static double Rot[3][3]={0.0};//
static double Rbo[3][3]={0.0};//
static double Rpo[3][3]={0.0};//模型上三点在osg下姿态
static double deltaPosition[3]={0.0};//


using namespace std;

 void GetDevicePosition(double P1[3],double P2[3])
 {
	 for(int i=0;i<3;i++)
	 {
		 P1[i]=0;//HumanStones.clsActiveDevices[0].vCenterPosition._v[i];
		 P2[i]=0;//HumanStones.clsActiveDevices[1].vCenterPosition._v[i];
	 }
 }
 void Product(const double * A1, const double * A2, double * A3)
 {
	 unsigned long  i,j,k;

	 for (i=0;i<3;i++) 
	 {
		 for(k=0;k<3;k++)
		 {
			 A3[i*3+k] = 0;
			 for (j=0;j<3;j++)
			 {
				 A3[i*3+k] += A1[i*3+j]*A2[j*3+k];
			 }
		 }
	 }
 }
 void EularToMatrix(const double Attitude[3],double  R[3][3])
 {
	 double S1,C1,S2,C2,S3,C3;
	 S1 = sin(Attitude[0]);
	 C1 = cos(Attitude[0]);
	 S2 = sin(Attitude[1]);
	 C2 = cos(Attitude[1]);
	 S3 = sin(Attitude[2]);
	 C3 = cos(Attitude[2]);

	 R[0][0] =  C2*C3 - S1*S2*S3;
	 R[0][1] =  C3*S1*S2 + C2*S3;
	 R[0][2] = -(C1*S2);
	 R[1][0] = -(C1*S3);
	 R[1][1] =  C1*C3;
	 R[1][2] =  S1;
	 R[2][0] =  C3*S2 + S1*C2*S3;
	 R[2][1] = -(C2*C3*S1)  + S2*S3;
	 R[2][2] =  C1*C2;

 }
 void EularExtract(const double Rbo[3][3],double Attitude[3])
 {
	 Attitude[0] =  asin(fLimitInOne(Rbo[1][2]));
	 if ( cos(Attitude[0]) > 0.0 )
		 Attitude[1] = atan2(-Rbo[0][2],Rbo[2][2]);
	 else
		Attitude[1]= atan2(Rbo[0][2],-Rbo[2][2]);
	 if ( cos(Attitude[0]) > 0.0 )
		Attitude[2]  = atan2(-Rbo[1][0],Rbo[1][1]);
	 else
		Attitude[2]  = atan2(Rbo[1][0],-Rbo[1][1]);
 }
 osg::ref_ptr<osg::PositionAttitudeTransform> rightpat; 
 osg::ref_ptr<osg::PositionAttitudeTransform> leftpat;
 double Ltb_tracker_right[3]={0.0};//tracker原点与骨头质心在tracker下
double Ltb_tracker_left[3]={0.0};//tracker原点与骨头质心在tracker下
 void RealTimeTracking(int Scale,const double Attitude[4],const double Position[3],const double Attitude2[4],const double Position2[3])
{
	double PositionInWorld[3];
	double AttitudeInWorld[3];
	double Position2InWorld[3];
	double Attitude2InWorld[3];
	double Rtc[3][3],Rbc[3][3];
	double Rct[3][3];
	struct SEularAngle  Angle;
	double Rtc0[3][3];

	double Q0[4]={0.0};
	double Rco0[3][3]={ { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
	double Rto0[3][3]={0.0};//tracker在osg坐标系下的姿态
	double Rpo0[3][3]={0.0};//骨头上标定点在osg坐标系姿态
	double Rbo0[3][3]={0.0};//骨头在osg坐标系姿态
	double relative_position_right[3];
	double relative_position_left[3];


	for (int i =0;i<4;i++)
		Q0[i] = Attitude[i];
	QuaternionToMatrix(Q0,Rtc0);//tracker在相机坐标系下的姿态
	double M1[3][3]={0.0};
	MatrixTranspose(3,3,&Rsp_right[0][0],&M1[0][0]);
	MatrixProduct(3,3,3,&M1[0][0],&Rtc0[0][0],&Rbo0[0][0]);
	//MatrixProduct(3,3,3,&Rtc0[0][0],&Rco0[0][0],&Rto0[0][0]);
	//MatrixProduct(3,3,3,&Rpt0[0][0],&Rto0[0][0],&Rpo0[0][0]);
	//MatrixProduct(3,3,3,&Rsp_right[0][0],&Rpo0[0][0],&Rbo0[0][0]);
	QuaternionExtract(Rbo0,Q0);
	Eular313Extract(Rbo0,&Angle);

	osg::Quat quat0;

	quat0._v[3]=Q0[0];
	quat0._v[0]=Q0[1];
	quat0._v[1]=Q0[2];
	quat0._v[2]=Q0[3];
	rightpat->setAttitude(quat0);
#if direct//直接用三轴转角
	osg::Quat quat0(Angle.Pitch,osg::Vec3d(0.0, 1.0, 0.0),Angle.Roll,osg::Vec3d(1.0, 0.0, 0.0),Angle.Yaw,osg::Vec3d(0.0, 0.0, 1.0));
	rightpat->setAttitude(quat0);
#endif	

	//加上相对位置
	MatrixTranspose(3,3,&Rtc0[0][0],&Rct[0][0]);
	MatrixProduct(3,3,1,&Rct[0][0],Ltb_tracker_right,relative_position_right);
	for (int i=0;i<3;i++)
	{
		relative_position_right[i]=0.0;
	}
	rightpat->setPosition(osg::Vec3d( Position[0]+relative_position_right[0],Position[1]+relative_position_right[1], Position[2]+relative_position_right[2]));


	double Q[4]={0.0};
	double Rco[3][3]={ { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
	double Rto[3][3]={0.0};//tracker在osg坐标系下的姿态
	double Rpo[3][3]={0.0};//骨头上标定点在osg坐标系姿态
	double Rbo[3][3]={0.0};//骨头在osg坐标系姿态
	//double Rsp[3][3]={0.0};

	for (int i =0;i<4;i++)
		Q[i] = Attitude2[i];
	QuaternionToMatrix(Q,Rtc);//tracker在相机坐标系下的姿态
	//QuaternionExtract(Rtc,Q);
	double M2[3][3]={0.0};
	MatrixTranspose(3,3,&Rsp_left[0][0],&M2[0][0]);
	MatrixProduct(3,3,3,&M2[0][0],&Rtc[0][0],&Rbo[0][0]);

	//MatrixProduct(3,3,3,&Rtc[0][0],&Rco[0][0],&Rto[0][0]);
	//MatrixProduct(3,3,3,&Rpt[0][0],&Rto[0][0],&Rpo[0][0]);
	//MatrixProduct(3,3,3,&Rsp_left[0][0],&Rpo[0][0],&Rbo[0][0]);

	QuaternionExtract(Rbo,Q);
	Eular313Extract(Rbo,&Angle);

	osg::Quat quat;
	quat._v[3]=Q[0];
	quat._v[0]=Q[1];
	quat._v[1]=Q[2];
	quat._v[2]=Q[3];
	leftpat->setAttitude(quat);
#if direct//直接用三轴转角
	osg::Quat quat1(Angle.Pitch,osg::Vec3d(0.0, 1.0, 0.0),Angle.Roll,osg::Vec3d(1.0, 0.0, 0.0),Angle.Yaw,osg::Vec3d(0.0, 0.0, 1.0));
	leftpat->setAttitude(quat1);
#endif

	//加上相对位置
	MatrixTranspose(3,3,&Rtc[0][0],&Rct[0][0]);
	MatrixProduct(3,3,1,&Rct[0][0],Ltb_tracker_left,relative_position_left);
	for (int i=0;i<3;i++)
	{
		relative_position_left[i]=0.0;
	}
	leftpat->setPosition(osg::Vec3d( Position2[0]+relative_position_left[0],Position2[1]+relative_position_left[1], Position2[2]+relative_position_left[2]));


 }
 /*char   buff[255]; */   
 const char* GetWorldIntersectPoint()
 {
	/* strcpy(buff,"100");*/
	 return worldCoordinate;
 }
void  GetTransformMatrix(const int id, const double P1[3],const double P2[3],const double Attitude[4],const double R1[3],const double R2[3],const double V1[3],const double V2[3],double M1[3][3],double M2[3][3])
 {
	 //P1为骨头实物上点在tracker下坐标系，P2为骨头模型上点在OSG坐标系下坐标，O为Tracker坐标系原点，O'为OSG坐标系原点
	//OP1为骨头实物上点到原点的向量（tracker坐标系），到投影到OSG坐标系下
	 //O'P2为骨头模型上点OSG原点的向量

	 //标定的时候能否得到Tracker姿态

	 //实物上某一点与tracker坐标系原点的位置矢量ox，得到ox在相机系下的投影ox'，对应的该点在osg下的位置o'x'，骨头质心在osg下的位置o'b
	//ob= ox'+x'b  tracker与骨头质心距离在osg/相机下的投影，将它投影到tracker坐标系下，为常值。

	//以后得到tracker姿态，得到相对位置在相机系投影，加上tracker的位置，认为是骨头质心位置。


	 double Rbi[3][3] = {0.0},Rei[3][3] = {0.0},Rop[3][3] = {0.0};

	 struct SEularAngle  Angle;
	 double Rtc[3][3];
	 QuaternionToMatrix(Attitude,Rtc);//tracker在相机坐标系下的姿态
	 if(id==0)
	 {
		 GetRdi(R1,R2,Rpt0);
		 GetRdi(V1,V2,Rpo0);
		 //模型在OSG下姿态
		double Q[4]={0.0};
		osg::Quat quat = rightpat->getAttitude();
		Q[0]=quat._v[3];//w();
		Q[1]=quat._v[0];//x();
		Q[2]=quat._v[1];//y();
		Q[3]=quat._v[2];//z();
		QuaternionToMatrix(Q,Rbo0);
		MatrixTranspose(3,3,&Rpo0[0][0],&Rop[0][0]);
		MatrixProduct(3,3,3,&Rbo0[0][0],&Rop[0][0],&Rsp_right[0][0]);
		 for(int i =0;i<3;i++)
		 {
			  for(int j =0;j<3;j++)
			  {
				  Rsp_right[i][j]=Rtc[i][j];
			  }
		 }

		//位置的处理 P1为实物选点，P2为osg对应的点
//		double Rtc[3][3]={0.0};//tracker在osg坐标系、相机下的姿态
		double Rct[3][3]={0.0};
		double Ltp_camera[3]={0.0};//tracker原点与p1位置在相机系下
		double bone_osg[3]={0.0};//骨头质心在osg下的位置
		double Lpb_osg[3]={0.0};//点与骨头质心在osg下的相对位置
		double Ltb_osg[3]={0.0};//tracker原点与骨头质心在osg下
		

		MatrixTranspose(3,3,&Rtc[0][0],&Rct[0][0]);
		MatrixProduct(3,3,1,&Rct[0][0],P1,Ltp_camera);
		osg::Vec3 position_osg = rightpat->getPosition();
		bone_osg[0]=position_osg._v[0];
		bone_osg[1]=position_osg._v[1];
		bone_osg[2]=position_osg._v[2];
		for (int i=0;i<3;i++)
		{
			Lpb_osg[i]= bone_osg[i]-P2[i];
		}
		for (int i=0;i<3;i++)
		{
			Ltb_osg[i]= Ltp_camera[i]+Lpb_osg[i];
		}
		MatrixProduct(3,3,1,&Rtc[0][0],Ltb_osg,Ltb_tracker_right);
	 }
	 if(id==1)
	 {
		 GetRdi(R1,R2,Rpt);
		 GetRdi(V1,V2,Rpo);
		//模型在OSG下姿态
		double Q[4]={0.0};
		osg::Quat quat = leftpat->getAttitude();
		Q[0]=quat._v[3];//w();
		Q[1]=quat._v[0];//x();
		Q[2]=quat._v[1];//y();
		Q[3]=quat._v[2];//z();
		QuaternionToMatrix(Q,Rbo);

		MatrixTranspose(3,3,&Rpo[0][0],&Rop[0][0]);
		MatrixProduct(3,3,3,&Rbo[0][0],&Rop[0][0],&Rsp_left[0][0]);

		 for(int i =0;i<3;i++)
		 {
			 for(int j =0;j<3;j++)
			 {
				 Rsp_left[i][j]=Rtc[i][j];
			 }
		 }


		//位置的处理 P1为实物选点，P2为osg对应的点
//		double Rtc[3][3]={0.0};//tracker在osg坐标系、相机下的姿态
		double Rct[3][3]={0.0};
		double Ltp_camera[3]={0.0};//tracker原点与p1位置在相机系下
		double bone_osg[3]={0.0};//骨头质心在osg下的位置
		double Lpb_osg[3]={0.0};//点与骨头质心在osg下的相对位置
		double Ltb_osg[3]={0.0};//tracker原点与骨头质心在osg下


		MatrixTranspose(3,3,&Rtc[0][0],&Rct[0][0]);
		MatrixProduct(3,3,1,&Rct[0][0],P1,Ltp_camera);
		osg::Vec3 position_osg = leftpat->getPosition();
		bone_osg[0]=position_osg._v[0];
		bone_osg[1]=position_osg._v[1];
		bone_osg[2]=position_osg._v[2];
		for (int i=0;i<3;i++)
		{
			Lpb_osg[i]= bone_osg[i]-P2[i];
		}
		for (int i=0;i<3;i++)
		{
			Ltb_osg[i]= Ltp_camera[i]+Lpb_osg[i];
		}
		MatrixProduct(3,3,1,&Rtc[0][0],Ltb_osg,Ltb_tracker_left);
	 }
	// MatrixProduct(3,3,3,&Rop[0][0],&Rpt[0][0],&Rot[0][0]);

	 for (int i=0;i<3;i++)
	 {
		 for (int j=0;j<3;j++)
		 {
			 M1[1][j]=Ltb_tracker_right[j];
		 }
	 }

	 for (int i=0;i<3;i++)
	 {
		 for (int j=0;j<3;j++)
		 {
			 M2[1][j]=Ltb_tracker_left[j];
		 }
	 }
 }
class RightRotateCallBack: public osg::NodeCallback
{
public:
	RightRotateCallBack():_rotateZ(0.0) {}

	virtual void operator()(osg::Node* node, osg::NodeVisitor* nv)
	{
		osg::PositionAttitudeTransform* pat = 
			dynamic_cast<osg::PositionAttitudeTransform*>(node);
		if(pat)
		{
			osg::Vec3 vec(0, 0, 1);

			//osg::Quat quat = osg::Quat(osg::DegreesToRadians(_rotateZ), osg::Z_AXIS);
			//pat->setAttitude(quat);
			//_rotatex=34;
			//osg::Quat quat2 = osg::Quat(osg::DegreesToRadians(_rotatex), osg::X_AXIS);	
			//pat->setAttitude(quat2*quat);
			//double position[3];
			//pat->setPosition(osg::Vec3d(_rotateZ*sin(_rotateZ/1000),0,0));
			//_rotateZ += 0.50;
		}

		traverse(node, nv);
	}

private:
	double _rotateZ;
	double _rotatex;
};
double rightangle;
class LeftRotateCallBack: public osg::NodeCallback
{
public:
	LeftRotateCallBack():_rotateZ(0.0) {}

	virtual void operator()(osg::Node* node, osg::NodeVisitor* nv)
	{
		osg::PositionAttitudeTransform* pat = 
			dynamic_cast<osg::PositionAttitudeTransform*>(node);
		if(pat)
		{
			osg::Vec3 vec(0, 0, 1);

			/*	osg::Quat quat = osg::Quat(osg::DegreesToRadians(_rotateZ), osg::Y_AXIS);
			pat->setAttitude(quat);
			_rotatex=34;
			osg::Quat quat2 = osg::Quat(osg::DegreesToRadians(_rotatex), osg::X_AXIS);	
			pat->setAttitude(quat2*quat);
			double position[3];
			pat->setPosition(osg::Vec3d(_rotateZ*sin(_rotateZ/1000),0,0));
			_rotateZ = 90;*/
		}

		traverse(node, nv);
	}

private:
	double _rotateZ;
	double _rotatex;
};
class InfoCallBack: public osg::NodeCallback
{
public:
	virtual void operator()(osg::Node* node, osg::NodeVisitor* nv)
	{
		osg::PositionAttitudeTransform* pat =
			dynamic_cast<osg::PositionAttitudeTransform*>(node);

		if(pat)
		{
			double angle = 0.0;
			osg::Vec3 axis;
			pat->getAttitude().getRotate(angle, axis);

			//std::cout << "Node is rotate around the axis(" << axis << "), "
			//	<<osg::RadiansToDegrees(angle) << "degrees" << std::endl;
		}

		traverse(node, nv);
	}
};

int osgRun( HWND m_hWnd,char * Path,char*Path2)

{
	RECT rect;
	::GetWindowRect(m_hWnd, &rect);
	root= new osg::Group; 
	viewer=new osgViewer::Viewer;	

	rightpat = new osg::PositionAttitudeTransform();
	leftpat = new osg::PositionAttitudeTransform();

	// 初始化图形环境GraphicsContext Traits
	osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits;
	osgViewer::GraphicsWindowWin32::WindowData* windata = new osgViewer::GraphicsWindowWin32::WindowData(m_hWnd);
	traits->x = 0;
	traits->y = 0;
	traits->width = rect.right - rect.left;
	traits->height = rect.bottom - rect.top;
	traits->windowDecoration = false;
	traits->doubleBuffer = true;
	traits->sharedContext = 0;
	traits->setInheritedWindowPixelFormat = true;
	traits->inheritedWindowData = windata;
	osg::ref_ptr<osg::GraphicsContext> gc = osg::GraphicsContext::createGraphicsContext(traits.get());
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
    root->addChild( g_grpMouse.get());

 

	osgText::Text* text = new osgText::Text; 
	root->addChild( createHUD_viewPoint( text));//加入HUD文字
	osg::ref_ptr< CHUD_viewPoint> pHUD= new CHUD_viewPoint( text);

	viewer->addEventHandler( pHUD.get());

	viewer->addEventHandler(new osgViewer::StatsHandler());
	viewer->addEventHandler(new osgViewer::WindowSizeHandler());


	osg::Node* rightmodel =  osgDB::readNodeFile(Path) ;//"cow.osg"


	rightpat->addChild(rightmodel);
	rightpat->setUpdateCallback(new RightRotateCallBack() );

	osg::Node* leftmodel =  osgDB::readNodeFile(Path2) ;//"cow.osg"

	leftpat->addChild(leftmodel);
	leftpat->setUpdateCallback(new LeftRotateCallBack() );


	root->addChild(rightpat.get());
	root->addChild(leftpat.get());

	viewer->addEventHandler(new PickHandler());
    viewer->setSceneData( root.get());

    viewer->realize();
    viewer->run() ;  
	
   // exit(0);
    return 0;

}
void EndDisplay()
{
/*	viewer->setSceneData(NULL);*/
	viewer->stopThreading();
	viewer->setDone(true);

	viewer->setQuitEventSetsDone(false);
	viewer->setKeyEventSetsDone(0x20);
}