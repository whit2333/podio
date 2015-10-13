// datamodel specific includes
#include "ExampleWithVectorMember.h"
#include "ExampleWithVectorMemberConst.h"
#include "ExampleWithVectorMemberObj.h"
#include "ExampleWithVectorMemberData.h"
#include "ExampleWithVectorMemberCollection.h"
#include <iostream>

ExampleWithVectorMember::ExampleWithVectorMember() : m_obj(new ExampleWithVectorMemberObj()){
 m_obj->acquire();
};



ExampleWithVectorMember::ExampleWithVectorMember(const ExampleWithVectorMember& other) : m_obj(other.m_obj) {
  m_obj->acquire();
}

ExampleWithVectorMember& ExampleWithVectorMember::operator=(const ExampleWithVectorMember& other) {
  if ( m_obj != nullptr) m_obj->release();
  m_obj = other.m_obj;
  return *this;
}

ExampleWithVectorMember::ExampleWithVectorMember(ExampleWithVectorMemberObj* obj) : m_obj(obj){
  if(m_obj != nullptr)
    m_obj->acquire();
}

ExampleWithVectorMember ExampleWithVectorMember::clone() const {
  return {new ExampleWithVectorMemberObj(*m_obj)};
}

ExampleWithVectorMember::~ExampleWithVectorMember(){
  if ( m_obj != nullptr) m_obj->release();
}

ExampleWithVectorMember::operator ConstExampleWithVectorMember() const {return ConstExampleWithVectorMember(m_obj);};



std::vector<int>::const_iterator ExampleWithVectorMember::count_begin() const {
  auto ret_value = m_obj->m_count->begin();
  std::advance(ret_value, m_obj->data.count_begin);
  return ret_value;
}

std::vector<int>::const_iterator ExampleWithVectorMember::count_end() const {
  auto ret_value = m_obj->m_count->begin();
  std::advance(ret_value, m_obj->data.count_end-1);
  return ++ret_value;
}

void ExampleWithVectorMember::addcount(int component) {
  m_obj->m_count->push_back(component);
  m_obj->data.count_end++;
}

bool  ExampleWithVectorMember::isAvailable() const {
  if (m_obj != nullptr) {
    return true;
  }
  return false;
}

const podio::ObjectID ExampleWithVectorMember::getObjectID() const {
  return m_obj->id;
}

bool ExampleWithVectorMember::operator==(const ConstExampleWithVectorMember& other) const {
     return (m_obj==other.m_obj);
}


//bool operator< (const ExampleWithVectorMember& p1, const ExampleWithVectorMember& p2 ) {
//  if( p1.m_containerID == p2.m_containerID ) {
//    return p1.m_index < p2.m_index;
//  } else {
//    return p1.m_containerID < p2.m_containerID;
//  }
//}
