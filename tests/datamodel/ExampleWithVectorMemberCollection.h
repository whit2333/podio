//AUTOMATICALLY GENERATED - DO NOT EDIT

#ifndef ExampleWithVectorMemberCollection_H
#define  ExampleWithVectorMemberCollection_H

#include <string>
#include <vector>
#include <deque>
#include <array>

// podio specific includes
#include "podio/ICollectionProvider.h"
#include "podio/CollectionBase.h"
#include "podio/CollectionIDTable.h"

// datamodel specific includes
#include "ExampleWithVectorMemberData.h"
#include "ExampleWithVectorMember.h"
#include "ExampleWithVectorMemberObj.h"

typedef std::vector<ExampleWithVectorMemberData> ExampleWithVectorMemberDataContainer;
typedef std::deque<ExampleWithVectorMemberObj*> ExampleWithVectorMemberObjPointerContainer;

class ExampleWithVectorMemberCollectionIterator {

  public:
    ExampleWithVectorMemberCollectionIterator(int index, const ExampleWithVectorMemberObjPointerContainer* collection) : m_index(index), m_object(nullptr), m_collection(collection) {}

    bool operator!=(const ExampleWithVectorMemberCollectionIterator& x) const {
      return m_index != x.m_index; //TODO: may not be complete
    }

    const ExampleWithVectorMember operator*() const;
    const ExampleWithVectorMember* operator->() const;
    const ExampleWithVectorMemberCollectionIterator& operator++() const;

  private:
    mutable int m_index;
    mutable ExampleWithVectorMember m_object;
    const ExampleWithVectorMemberObjPointerContainer* m_collection;
};

/**
A Collection is identified by an ID.
*/

class ExampleWithVectorMemberCollection : public podio::CollectionBase {

public:
  typedef const ExampleWithVectorMemberCollectionIterator const_iterator;

  ExampleWithVectorMemberCollection();
//  ExampleWithVectorMemberCollection(const ExampleWithVectorMemberCollection& ) = delete; // deletion doesn't work w/ ROOT IO ! :-(
//  ExampleWithVectorMemberCollection(ExampleWithVectorMemberVector* data, int collectionID);
  ~ExampleWithVectorMemberCollection(){};

  void clear();
  /// Append a new object to the collection, and return this object.
  ExampleWithVectorMember create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given
  template<typename... Args>
  ExampleWithVectorMember create(Args&&... args);
  int size() const;

  /// Returns the object of given index
  const ExampleWithVectorMember operator[](int index) const;

  /// Append object to the collection
  void push_back(ConstExampleWithVectorMember object);

  void prepareForWrite();
  void prepareAfterRead();
  void setBuffer(void* address);
  bool setReferences(const podio::ICollectionProvider* collectionProvider);

  podio::CollRefCollection* referenceCollections() { return m_refCollections;};

  void setID(unsigned ID){m_collectionID = ID;};

  // support for the iterator protocol
  const const_iterator begin() const {
    return const_iterator(0, &m_entries);
  }
  const	const_iterator end() const {
    return const_iterator(m_entries.size(), &m_entries);
  }

  /// returns the address of the pointer to the data buffer
  void* getBufferAddress() { return (void*)&m_data;};

  /// returns the pointer to the data buffer
  std::vector<ExampleWithVectorMemberData>* _getBuffer() { return m_data;};

   

private:
  int m_collectionID;
  ExampleWithVectorMemberObjPointerContainer m_entries;
  // members to handle 1-to-N-relations

  // members to handle streaming
  podio::CollRefCollection* m_refCollections;
  ExampleWithVectorMemberDataContainer* m_data;
};

template<typename... Args>
ExampleWithVectorMember  ExampleWithVectorMemberCollection::create(Args&&... args){
  int size = m_entries.size();
  auto obj = new ExampleWithVectorMemberObj({size,m_collectionID},{args...});
  m_entries.push_back(obj);
  return ExampleWithVectorMember(obj);
}



#endif
