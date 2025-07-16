#ifndef LSMS_LUASUPPORT_H
#define LSMS_LUASUPPORT_H

#include "Real.hpp"
#include "lua.hpp"

void luaStackDump(lua_State *L);

bool luaGetStrN(lua_State *L, const char *name, char *s, int n);
bool luaGetReal(lua_State *L, const char *name, Real *val);
bool luaGetInteger(lua_State *L, const char *name, int *val);
bool luaGetInteger(lua_State *L, const char *name, unsigned int *val);
bool luaGetBoolean(lua_State *L, const char *name, bool *val);
bool luaGetFieldInTable(lua_State *L, const char *name, const char *field);
bool luaGetFieldFromStack(lua_State *L, const char *field);
bool luaGetPositionInTable(lua_State *L, const char *name, int idx);
bool luaGetRealFieldInTable(lua_State *L, const char *name, const char *field,
                            Real *val);
bool luaGetIntegerFieldInTable(lua_State *L, const char *name,
                               const char *field, int *val);
bool luaGetRealPositionInTable(lua_State *L, const char *name, int idx,
                               Real *val);
bool luaGetIntegerPositionInTable(lua_State *L, const char *name, int idx,
                                  int *val);
bool luaGetStrNFromStack(lua_State *L, const char *name, char *s, int n);

template <typename T>
bool luaGetRealPositionFromStack(lua_State *L, int idx, T *val) {
  lua_pushinteger(L, idx);
  lua_gettable(L, -2);
  if (!lua_isnumber(L, -1)) {
    lua_pop(L, 1);
    return false;
  }
  *val = lua_tonumber(L, -1);
  lua_pop(L, 1);
  return true;
}


bool luaGetIntegerPositionFromStack(lua_State *L, int idx, int *val);
bool luaGetRealFieldFromStack(lua_State *L, const char *field, Real *val);
bool luaGetIntegerFieldFromStack(lua_State *L, const char *field, int *val);

#endif
