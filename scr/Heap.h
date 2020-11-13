#pragma once
#define DefaltSize 1000

template <class Type>
class Heap{
public:
	virtual bool Insert(const Type& )=0;
	virtual Type * Delete(Type& )=0;
};
// MinHeap the 
template <class Type>
class Min_Heap:public Heap<Type>
{
public:
	bool IsFull();// if full return true
	bool IsEmpty();// if empty return true
	bool EnLarge();
	bool Contract();
	int GetNum() {return n;}
	int GetMaxSize() {return MaxSize;}
	bool Insert(const Type& x);
	Type* Delete(Type& x);
	Type Pop();
    Min_Heap(int sz=DefaltSize);
	Min_Heap(const Min_Heap& other);
	Min_Heap& operator=(const Min_Heap& other);
	~Min_Heap() { delete[] storge; MaxSize=n=0; storge=NULL;};
//private:
public:
	Type* storge;
	int n;//the number of the elements in the min heap 
	int MaxSize;//the upperbound of the min heap
};

template <class Type>
Min_Heap<Type>::Min_Heap(int sz)
{
	 n=0;
     MaxSize=sz;
	 storge=new Type[MaxSize+1];
	 memset(storge,0,(MaxSize+1)*sizeof(Type));
}
template <class Type>
Min_Heap<Type>::Min_Heap(const Min_Heap<Type> &other)
{
	this->n=other.n;
	this->MaxSize=other.MaxSize;
	this->storge=new Type[MaxSize+1];
	memcpy(storge,other.storge,sizeof(Type)*(MaxSize+1));
}
template <class Type>
Min_Heap<Type>& Min_Heap<Type>::operator =(const Min_Heap<Type>& other)
{
	//检查自赋值
	if(this==&other)
	{
		return (*this);
	}
    //清除原有的内存
	delete[] storge;
	this->n=other.GetNum();
	this->MaxSize=other.GetMaxSize();
	storge=new Type[MaxSize+1];
	memcpy(storge,other.storge,sizeof(Type)*(MaxSize+1));
	return (*this);
}
template <class Type>
bool Min_Heap<Type>::IsFull()
{
	return n==MaxSize;
}

template <class Type>
bool Min_Heap<Type>::IsEmpty()
{
	return n==0;
}

template <class Type>
bool Min_Heap<Type>::EnLarge()
{
    MaxSize*=2;
	Type * storge_temp=new Type[MaxSize+1];
	memcpy(storge_temp,storge,(n+1)*sizeof(Type));
	delete[] storge;
	storge=storge_temp;
	storge_temp=NULL;
	return true;
}
template <class Type>
bool Min_Heap<Type>::Contract()
{
	if (n>(MaxSize))
	{
		return false;
	}
	MaxSize/=2;
	Type* storge_temp=new Type[MaxSize+1];
	memcpy(storge_temp,storge,(n+1)*sizeof(Type));
	delete[] storge;
	storge=storge_temp;
	storge_temp=NULL;
	return true;
}
template <class Type>
bool Min_Heap<Type>::Insert(const Type &x)
{
	if (this->IsFull())
	{
		//return false;
		this->EnLarge();
	}
	n++;
	int i=n;
	for(;i>1;)
	{
       if (x >= storge[i/2])
       {
		   break;
       }
	   storge[i]=storge[i/2];
	   i=i/2;
	}
	storge[i]=x;
	return true;
}

template <class Type>
Type* Min_Heap<Type>::Delete(Type& x)
{
   if (this->IsEmpty())
   {
	   return NULL;
   }
   x=this->storge[1];
   Type y=storge[n--];
   int i=1,j=2;
   for (;j<=n;)
   {
	   if (j<n)
	   {
		   if (storge[j+1]<storge[j])
		   {
			   j++;
		   }
	   }
	   if (y<storge[j])
	   {
		   break;
	   }
	   storge[i]=storge[j];
	   i=j;
	   j*=2;
   }
   storge[i]=y;
   if (n<(MaxSize/4))
   {
	   this->Contract();
   }
   return &x;
}
template<class Type>
Type Min_Heap<Type>::Pop()
{
	Type x;
 	if(this->IsEmpty())
	{
       return x;
	}
	x=this->storge[1];
	Type y=storge[n--];
	int i=1,j=2;
	for (;j<=n;)
	{
		if (j<n)
		{
			if(storge[j+1]<storge[j])
			{
				j++;
			}
		}
		if(y<storge[j])
		{
			break;
		}
		storge[i]=storge[j];
		i=j;
		j*=2;
	}
	storge[i]=y;
	return x;
}
//Min_Heap_Pointer: storge the pointer of the element
 template<class Type>
 class Min_Heap_Pointer:public Heap<Type>
 {
 public:
	bool IsFull() {return n==MaxSize;}
	bool IsEmpty() {return n==0;}
 	bool EnLarge();
 	bool Insert(const Type& x);
	Type* Delete(Type& x);
	Type Pop();
	Min_Heap_Pointer(int sz=DefaltSize);
	Min_Heap_Pointer(const Min_Heap_Pointer& other);
	Min_Heap_Pointer& operator=(const Min_Heap_Pointer& other);
	~Min_Heap_Pointer() {delete[] storge;n=0;MaxSize=0;}
 //private:
 public:
	 Type* storge;
	 int n;
	 int MaxSize;
 };
 template <class Type>
 Min_Heap_Pointer<Type>::Min_Heap_Pointer(int sz)
 {
	 this->n=0;
	 this->MaxSize=sz;
	 storge=new Type[MaxSize+1];
	 memset(storge,0,(MaxSize+1)*sizeof(Type));
 }
template <class Type>
bool Min_Heap_Pointer<Type>::EnLarge()
{
	MaxSize*=2;
	Type* storge_temp=new Type[MaxSize+1];
	memcpy(storge_temp,storge,(n+1)*sizeof(Type));
	delete[] storge;
	storge=storge_temp;
	storge_temp=NULL;
	return true;
}
template <class Type>
bool Min_Heap_Pointer<Type>::Insert(const Type& x)
{
	if(this->IsFull())
	{
		return false;
	}
	n++;
	int i=n;
	for(;i>1;)
	{
		if ((*x)>=(*storge[i/2]))//pointer
		{
			break;
		}
		storge[i]=storge[i/2];
		i/=2;
	}
    storge[i]=x;
	return true;
}
template <class Type>
Type* Min_Heap_Pointer<Type>::Delete(Type &x)
{
	if(this->IsEmpty())
	{
		return NULL;
	}
	x=storge[1];
	Type y=storge[n--];
    int i=1,j=2;
	for (;j<=n;)
	{
		if (j<n)
		{
			if((*storge[j+1])<(*storge[j]))
			{
				j++;
			}
		}
		if ((*y)<(*storge[j]))
		{
			break;
		}
		storge[i]=storge[j];
		i=j;
		j*=2;
	}
	storge[i]=y;
	return &x;
}
template <class Type>
Type Min_Heap_Pointer<Type>::Pop()
{
    Type x;
	if (this->IsEmpty())
	{
		return x;
	}
	x=storge[1];
	Type y=storge[n--];
	int i=1,j=2;
	for (;j<=n;)
	{
		if (j<n)
		{
			if ((*storge[j+1])<(*storge[j]))
			{
				j++;
			}
		}
		if ((*y)<(*storge[j]))
		{
			break;
		}
        storge[i]=storge[j];
		i=j;
		j*=2;
	}
	storge[i]=y;
	return x;
}

template <class Type>
Min_Heap_Pointer<Type>::Min_Heap_Pointer(const Min_Heap_Pointer<Type> &other)
{
	this->n=other.n;
	this->MaxSize=other.MaxSize;
	this->storge=new Type[MaxSize+1];
	memcpy(storge,other.storge,(n+1)*sizeof(Type));
}
template <class Type>
Min_Heap_Pointer<Type>& Min_Heap_Pointer<Type>::operator =(const Min_Heap_Pointer<Type> &other)
{
	//检查自赋值
	if(this=&other)
	{
		return (*this);
	}
	this->n=other.n;
	this->MaxSize=other.MaxSize;
	delete[] storge;
	storge=new Type[MaxSize+1];
	memcpy(storge,other.storge,sizeof(Type)*(n+1));
	return (*this);
}