#pragma once
#include "vector"
#include <memory>
#include "iostream"

class LinearMemResource {
public:
    LinearMemResource(): cur_pos(0), cur_size(0) {
    }

    ~LinearMemResource() {
        for (auto w: buffers) {
            ::operator delete(w);
        }
    }

    void* allocate(std::size_t size) {
        if (size + cur_pos > cur_size) {
            add_new_buffer();
        }
        auto res = reinterpret_cast<char*>(buffers.back()) + cur_pos;
        cur_pos += size;
        return res;
    }

    void clear() {
        if (buffers.size() < 2) {
            cur_pos = 0;
            return;
        }
        auto last = buffers.back();
        for (int i = 0; i < static_cast<int>(buffers.size()) - 1; ++i) {
            ::operator delete(buffers[i]);
        }
        buffers.clear();
        buffers.push_back(last);
        cur_pos = 0;
    }
private:
    void add_new_buffer() {
        if (buffers.empty()) {
            cur_size = (1 << 15);
        } else {
            cur_size *= 2;
        }
        buffers.push_back(::operator new(cur_size));
        cur_pos = 0;
        buffers_added += 1;
//        std::cout << "buffers added: " << buffers_added << std::endl;
    }

    int buffers_added=0;
    std::size_t cur_pos;
    std::size_t cur_size;
    std::vector<void*> buffers;
};

template<class T>
class StackAllocator {
    template<class U>
    friend class StackAllocator;
public:
    using value_type = T;

    StackAllocator(): resource_(std::make_shared<LinearMemResource>()) {
    }

    template<class U>
    StackAllocator(const StackAllocator<U>& other): resource_(other.resource_) {
    }

    T* allocate(std::size_t n) {
        return (T*) resource_->allocate(sizeof(T) * n);
    }

    void deallocate(T*, std::size_t) {
    }

    void consolidate() {
        resource_->clear();
    }

private:
    std::shared_ptr<LinearMemResource> resource_;
};