#pragma once
#include <assert.h>

template<typename T>
class Singleton
{
public:
	static inline void Create()
	{
		if(!pInstance)
			new T();
	}
	static inline T* Get()
	{
		assert(pInstance);
		return pInstance;
	}
	
protected:

	Singleton()
	{
		assert(!pInstance);
		pInstance = static_cast<T*>(this);
	}

	Singleton& operator= (const Singleton&);

private:
	static	T* pInstance;

};

template <typename T> T* Singleton <T>::pInstance = 0;
