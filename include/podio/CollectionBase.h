#ifndef COLLECTIONBASE_H
#define COLLECTIONBASE_H

#include <string>
#include <utility>
#include <vector>

#include "TTree.h"

#include "podio/ObjectID.h"

namespace podio {
  // forward declarations
  class ObjectID;
  class ICollectionProvider;
  class CollectionBase;

  typedef std::vector<std::pair<std::string,podio::CollectionBase*>> CollRegistry;
  typedef std::vector<std::vector<podio::ObjectID>*> CollRefCollection;

  template<class T>
  class UnderlyingTypeStuff {
  public:
    typedef T object_type;
    T* getBufferAddress2(){
      return static_cast<T*>(this->getBufferAddress());
    }
  };


  //class CollectionBuffer {
  //public:
  //  void* data;
  //  CollRefCollection* references;
  //};

  class CollectionBase {
  public:

    virtual void makeBranch(TTree* t, std::string name, std::string collClassName) = 0;
      //m_datatree->Branch(name.c_str(),  collClassName.c_str(), coll->getBufferAddress(), 32000, 199);

    /// prepare buffers for serialization
    virtual void  prepareForWrite() = 0;

    //virtual void  write(CollectionBuffer& buffer) = 0;
    //virtual void  read(CollectionBuffer& buffer) = 0;
    
    /// re-create collection from buffers after read
    virtual void  prepareAfterRead() = 0;

    /// initialize references after read
    virtual bool  setReferences(const ICollectionProvider* collectionProvider) = 0;

    /// set collection ID
    virtual void  setID(unsigned id) = 0;

    /// set I/O buffer
    virtual void  setBuffer(void*) = 0;

    /// get address of the pointer to the I/O buffer
    virtual void* getBufferAddress() = 0;

    /// check for validity of the container after read
    virtual bool isValid() const = 0;

    /// destructor
    virtual ~CollectionBase(){};

    /// clear the collection and all internal states
    virtual void clear() = 0 ;

    /// return the buffers containing the object-relation information
    virtual CollRefCollection* referenceCollections() = 0;
  };

} // namespace

#endif
