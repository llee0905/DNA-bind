/*****************************************************************************/
/***************** Copyright (C) 2020-2021, Richard Spinney. *****************/
/*****************************************************************************/
//                                                                           //
//    This program is free software: you can redistribute it and/or modify   //
//    it under the terms of the GNU General Public License as published by   //
//    the Free Software Foundation, either version 3 of the License, or      //
//    (at your option) any later version.                                    //
//                                                                           //
//    This program is distributed in the hope that it will be useful,        //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of         //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          //
//    GNU General Public License for more details.                           //
//                                                                           //
//    You should have received a copy of the GNU General Public License      //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef PYTHON_INTERFACE_H
#define PYTHON_INTERFACE_H

#include <string>
#include <vector>
#include <stdexcept>
#include <any>
#include <Python.h>

//interface to the python c-API for calling python functions and getting the return result using std::any functionality - requires c++17

//NOTE: only supports simple data types (int,double,char*), [not references or pointers and not const] & std::strings, and vectors of them
//I.e. can't currently do maps, lists etc.

class PythonInterface{
private:
    PyGILState_STATE m_gilState;
    PyObject *m_pName,*m_pModule,*m_pFunc,*m_pArgs,*m_pValue;    
    std::vector<PyObject*> m_argVec;

    template<typename T> PyObject* convertToPython(const std::vector<T> &a_vec) const {return vectorToList(a_vec);}
    PyObject* convertToPython(const char* a_cStr) const {return PyUnicode_FromString(a_cStr);}
    PyObject* convertToPython(const std::string &a_str) const {return PyUnicode_FromString(a_str.c_str());}
    template<typename T> PyObject* convertToPython(const T &a_type) const{
        static_assert(std::is_floating_point<T>::value || std::is_integral<T>::value,"Type not supported for convert_to_python function");
        if constexpr(std::is_floating_point<T>::value) return PyFloat_FromDouble(a_type);
        else if constexpr(std::is_integral<T>::value) return PyLong_FromLong(a_type);
        else throw std::logic_error("convert to python: unsupported type"); //shouldn't be possible to get here.        
    }

    std::any convertFromPython(PyObject* a_obj) const{
        if (PyLong_Check(a_obj)){
            return PyLong_AsLong(a_obj);    // c++17 converts it to a std::any automatically - same for all below
        }
        else if (PyFloat_Check(a_obj)){
            return PyFloat_AsDouble(a_obj);
        }
        else if (PyUnicode_Check(a_obj)){
            return std::string(PyUnicode_AsUTF8AndSize(a_obj,NULL)); //2nd arg: Py_ssize_t* size = NULL;
        }
        else if(PyList_Check(a_obj)){
            std::vector<std::any> data;
            for(Py_ssize_t listIndex = 0; listIndex < PyList_Size(a_obj); ++listIndex) {
                PyObject *value = PyList_GetItem(a_obj,listIndex);
                data.push_back(convertFromPython(value)); //note the recursion
            }
            return data;
        }
        else if(PyTuple_Check(a_obj)){
            std::vector<std::any> data;
            for(Py_ssize_t tupleIndex = 0; tupleIndex < PyTuple_Size(a_obj); ++tupleIndex) {
                PyObject *value = PyTuple_GetItem(a_obj,tupleIndex);
                data.push_back(convertFromPython(value));
            }
            return data;
        }
        else{
            throw std::logic_error("convert from python: unsupported type\n");     
        }
    }    
    
    template<typename T> 
    PyObject* vectorToList(const std::vector<T> &a_vec) const{
        PyObject* listObj = PyList_New(a_vec.size());
        if (!listObj) throw std::logic_error("Unable to allocate memory for Python list\n");
        for (uint32_t vecIndex = 0; vecIndex < a_vec.size(); ++vecIndex){
            PyObject *val = convertToPython(a_vec[vecIndex]);
            if (!val) {
                Py_DECREF(listObj);
                throw std::logic_error("Unable to allocate memory for Python list\n");
            }
            PyList_SET_ITEM(listObj,vecIndex,val); // "steals" the reference val, no need to Py_DECREF() it.
                                                   // note alternative func: PyTuple_SET_ITEM for tuples
        }
        return listObj;
    }
    
    //recursive variadic template pattern for processing variable number of arguments
    template <typename T>
    void packageArguments(const T &a_type){ // last item in list 
        m_argVec.push_back(convertToPython(a_type));
        if (!m_argVec.back()) {
            Py_DECREF(m_pModule);
            for (auto arg : m_argVec) Py_DECREF(arg);
            throw std::logic_error("cannot convert argument\n");
        }
    }
    
    template <typename T,typename... Args>
    void packageArguments(const T &a_type,const Args&... a_args){
        m_argVec.push_back(convertToPython(a_type)); // deal with first in list
        if (!m_argVec.back()) {            
            Py_DECREF(m_pModule);
            for (auto arg : m_argVec) Py_DECREF(arg);
            throw std::logic_error("cannot convert argument\n");
        }
        packageArguments(a_args...); //send the rest to get processed
    }

public:

    //public interface: pass name of module that contains function, name of function, then list of arguments
    // returns a std::any type equal to the return type of the computation where all lists/tuples are turned into vectors
    // the caller is responsible for std::any_cast<> -ing the result into the desired c++ type

    template<typename... Args>
    std::any callFunction(const std::string &a_module,const std::string &a_func,const Args&... a_args){
    	m_argVec.clear(); //just in case
        m_gilState = PyGILState_Ensure();
        m_pName = PyUnicode_DecodeFSDefault(a_module.c_str());// modeule/script where the function is e.g. nupack_functions
        m_pModule = PyImport_Import(m_pName);
        Py_DECREF(m_pName);
        if (m_pModule != NULL) {
            m_pFunc = PyObject_GetAttrString(m_pModule,a_func.c_str());// the function to be called
            if (m_pFunc && PyCallable_Check(m_pFunc)) {
                packageArguments(a_args...); //fill placeholder vector of python objects equal to the passed parameters args...  in order.
                m_pArgs = PyTuple_New(m_argVec.size()); // create a python tuple w/ size equal to the now known number of parameters
                for (size_t argIndex=0; argIndex<m_argVec.size(); ++argIndex) //fill the tuple with the arguments
                    PyTuple_SetItem(m_pArgs, argIndex, m_argVec[argIndex]); //"steals" the reference m_argVec[i], no need to Py_DECREF() it.
                m_argVec.clear(); //clear the vector of arguments
                m_pValue = PyObject_CallObject(m_pFunc,m_pArgs); //call the function
                Py_DECREF(m_pArgs); //release the arguments
                if (!m_pValue) {
                    Py_DECREF(m_pFunc);
                    Py_DECREF(m_pModule);
                    PyErr_Print();
                    PyGILState_Release(m_gilState);
                    throw std::logic_error("Call failed\n");
                }
            }
            else {
                if (PyErr_Occurred()) PyErr_Print();
                Py_DECREF(m_pFunc);
                Py_DECREF(m_pModule);
                PyGILState_Release(m_gilState);
                throw std::logic_error("Cannot find function\n");
            }
            Py_XDECREF(m_pFunc); // release the function should it still exist
            Py_DECREF(m_pModule);  // release the module/script
        }
        else {
            PyErr_Print();
            PyGILState_Release(m_gilState);
            throw std::logic_error("Failed to find script/module\n");
        }
        std::any returnVal = convertFromPython(m_pValue); //pass the return value of the function back as a std::any
        Py_DECREF(m_pValue); //release the python return object       
        PyGILState_Release(m_gilState);
        return returnVal; 
        // up to the caller to translate the returned std::any into the structures they require
        // i.e. they need to know what the python function returns
        // note: *everything* is a std::any
        // i.e. if the return type is a vector it will first be need to be any_cast to std::vector<std::any>, 
        // and then all elements will also need to be any_cast too
    }   
    
    PythonInterface(const std::string a_pathToPython){
        setenv("PYTHONPATH",a_pathToPython.c_str(),1);
        Py_Initialize();
    }
    
    ~PythonInterface(){
        //Py_FinalizeEx(); // finalize python interpreter     
                           // bug in NumPy where it doesn't release and so crashes on next Py_Initialise()
                           // see, for example a) https://github.com/numpy/numpy/issues/8097
                           //                  b) https://github.com/numpy/numpy/issues/11925
                           // "solution" is to just not call Py_Finalise()
    }
};

#endif /*PYTHON_INTERFACE_H*/
