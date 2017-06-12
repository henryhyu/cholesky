// empty, C/C++ too hard

#ifndef _SOLVER_H
#define _SOLVER_H

void put_json_entry(json_object *jobj, const char *name, const char *ylab, const char *xlab);
void update_ydata(json_object *jobj, double value, const char *name, const char *fname)
void update_xdata(json_object *jobj, double value, const char *name)
#endif