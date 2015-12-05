#include <stdio.h>
#define DebugPrintf printf

#define DebugAlltitude 0x0001
#define DebugOrbit 0x0002
#define DebugControlMode 0x04
#define DebugACSStatus 0x08
#define DebugCraft 0x10
#define DebugACS 0x20
#define DebugFault 0x40
#define DebugTime 0x80

#define DebugALL 0xffff
//#define DebugMode (DebugTime|DebugFault)
#define DebugMode (DebugAlltitude | DebugControlMode | DebugACSStatus |DebugCraft)